//! `zipper`: Merge unmapped and mapped BAM files with metadata transfer
//!
//! This module implements the `zipper` tool, which merges information from an unmapped BAM
//! (typically from FASTQ conversion) with a mapped BAM (after alignment), transferring all
//! tags and metadata from the unmapped reads to the mapped reads.
//!
//! # Requirements
//!
//! Both input BAMs MUST be:
//! - Queryname sorted or grouped (all records with same name are consecutive)
//! - Have the same queryname ordering
//!
//! # Workflow
//!
//! 1. Read templates (groups of records with same name) from both BAMs in parallel
//! 2. For each template:
//!    - Remove specified tags from mapped reads
//!    - Copy all tags from unmapped to mapped reads
//!    - Apply reverse/revcomp transformations for negative strand reads
//!    - Transfer QC pass/fail flags
//! 3. Write merged records to output
//!
//! # Tag Manipulation
//!
//! - **Remove**: Tags listed in `--tags-to-remove` are removed from mapped reads
//! - **Reverse**: Tags listed in `--tags-to-reverse` are reversed for negative strand reads
//! - **Revcomp**: Tags listed in `--tags-to-revcomp` are reverse complemented for negative strand
//!
//! Named tag sets (e.g., "Consensus") are expanded to their constituent tags.
//!
//! # Example
//!
//! ```bash
//! fgumi zipper \
//!   -i mapped.bam \
//!   -u unmapped.bam \
//!   -r reference.fa \
//!   -o output.bam \
//!   --tags-to-reverse Consensus \
//!   --tags-to-revcomp Consensus
//! ```
//!
//! # Architecture Notes
//!
//! - `Template`: Represents a group of records with the same read name
//! - `TemplateIterator`: Lazily groups records by name from a BAM reader
//! - `TagInfo`: Holds sets of tags to remove/reverse/revcomp
//! - `merge()`: Core function that transfers metadata between templates
use crate::commands::command::Command;
use crate::commands::common::CompressionOptions;
use anyhow::{Context, Result};
use clap::Parser;
use fgumi_lib::bam_io::{create_bam_reader, create_bam_writer};
use fgumi_lib::batched_sam_reader::BatchedSamReader;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::sam::{
    buf_value_to_smallest_signed_int, check_sort, revcomp_buf_value, reverse_buf_value,
};
use fgumi_lib::template::{Template, TemplateIterator};
use fgumi_lib::umi::TagInfo;
use fgumi_lib::validation::validate_file_exists;
use log::{debug, info};
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::record_buf::RecordBuf;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Command-line arguments for `zipper`
///
/// Merges unmapped and mapped BAM files, transferring tags and metadata from unmapped
/// reads to their corresponding mapped reads. Supports tag manipulation (removal, reversal,
/// reverse complement) based on read strand.
#[derive(Parser, Debug)]
#[command(
    name = "zipper",
    author,
    version,
    about = "\x1b[38;5;72m[ALIGNMENT]\x1b[0m      \x1b[36mZip unmapped BAM with aligned BAM\x1b[0m",
    long_about = r#"
Merges unmapped and mapped BAM files, transferring tags and metadata.

Takes an unmapped BAM (typically from FASTQ) and a mapped BAM (after alignment) and merges
them, copying tags from the unmapped to mapped reads. Both BAMs must be queryname sorted or
grouped, and have the same read name ordering.

The tool transfers tags from the unmapped reads to their corresponding mapped reads. For reads
mapped to the negative strand, tags can be optionally reversed or reverse-complemented. All
QC pass/fail flags are also transferred from the unmapped to mapped reads.

## Tag Manipulation

You can specify which tags to manipulate for reads mapped to the negative strand:
- --tags-to-reverse: Reverses array and string tags (e.g., [1,2,3] becomes [3,2,1])
- --tags-to-revcomp: Reverse complements sequence tags (e.g., AGAGG becomes CCTCT)

Named tag sets like "Consensus" are automatically expanded to their constituent tags:
- Consensus: aD bD cD aM bM cM aE bE cE ad bd cd ae be ce ac bc

## Default Behavior

By default, input is read from stdin and output is written to stdout, allowing for streaming
workflows like:

  fgumi zipper -u unmapped.bam -r reference.fa | samtools sort -@ $(nproc) -o output.bam
"#
)]
#[command(verbatim_doc_comment)]
pub struct Zipper {
    /// Input mapped SAM file (or `-` for stdin)
    #[arg(short = 'i', long, default_value = "-")]
    pub input: PathBuf,

    /// Input unmapped BAM file containing original tags
    #[arg(short = 'u', long)]
    pub unmapped: PathBuf,

    /// Reference FASTA file (must have accompanying .dict file)
    #[arg(short = 'r', long)]
    pub reference: PathBuf,

    /// Output BAM file (or `-` for stdout)
    #[arg(short = 'o', long, default_value = "-")]
    pub output: PathBuf,

    /// Tags to remove from mapped reads before copying unmapped tags
    #[arg(long, value_delimiter = ',')]
    pub tags_to_remove: Vec<String>,

    /// Tags to reverse for reads mapped to negative strand
    #[arg(long, value_delimiter = ',')]
    pub tags_to_reverse: Vec<String>,

    /// Tags to reverse complement for reads mapped to negative strand
    #[arg(long, value_delimiter = ',')]
    pub tags_to_revcomp: Vec<String>,

    /// Buffer size for template channel (default: 50000)
    #[arg(short = 'b', long, default_value = "50000")]
    pub buffer: usize,

    /// Number of threads to use for processing (default: 1, single-threaded)
    #[arg(long, short = 't', default_value = "1")]
    pub threads: usize,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// BWA -K parameter value (bases per batch). Used to optimize buffer sizing
    /// for stdin input. The buffer grows adaptively based on observed bytes per batch.
    /// Default matches common bwa mem usage.
    #[arg(short = 'K', long = "bwa-chunk-size", default_value = "150000000")]
    pub bwa_chunk_size: u64,

    /// Exclude reads from the unmapped BAM that are not present in the aligned BAM.
    /// Useful when reads were intentionally removed (e.g., by adapter trimming) prior to alignment.
    #[arg(long = "exclude-missing-reads", default_value = "false")]
    pub exclude_missing_reads: bool,
}

/// Builds the output BAM header from unmapped and mapped headers
///
/// Merges information from both input headers:
/// - Reference sequences from the reference dictionary
/// - Comments from both input files
/// - Read groups from both inputs (mapped takes precedence)
/// - Program records from both inputs (mapped takes precedence)
/// - Header lines from unmapped BAM
///
/// # Arguments
///
/// * `unmapped` - Header from the unmapped BAM
/// * `mapped` - Header from the mapped BAM  
/// * `dict_path` - Path to reference FASTA (will open .dict file)
///
/// # Returns
///
/// Combined header for the output BAM
///
/// # Errors
///
/// Returns an error if the reference dictionary cannot be opened or parsed
pub fn build_output_header(unmapped: &Header, mapped: &Header, dict_path: &Path) -> Result<Header> {
    let dict_file = std::fs::File::open(dict_path.with_extension("dict"))
        .context("Failed to open reference dictionary")?;
    let mut dict_reader = noodles::sam::io::Reader::new(std::io::BufReader::new(dict_file));
    let dict_header = dict_reader.read_header()?;

    let mut header = Header::builder();

    for (name, map) in dict_header.reference_sequences() {
        header = header.add_reference_sequence(name.clone(), map.clone());
    }

    for comment in unmapped.comments() {
        header = header.add_comment(comment.clone());
    }
    for comment in mapped.comments() {
        header = header.add_comment(comment.clone());
    }

    let mut rg_ids = HashSet::new();
    for (id, rg) in mapped.read_groups() {
        header = header.add_read_group(id.clone(), rg.clone());
        rg_ids.insert(id.clone());
    }
    for (id, rg) in unmapped.read_groups() {
        if !rg_ids.contains(id) {
            header = header.add_read_group(id.clone(), rg.clone());
        }
    }

    let mut pg_ids = HashSet::new();
    for (id, pg) in mapped.programs().as_ref() {
        header = header.add_program(id.clone(), pg.clone());
        pg_ids.insert(id.clone());
    }
    for (id, pg) in unmapped.programs().as_ref() {
        if !pg_ids.contains(id) {
            header = header.add_program(id.clone(), pg.clone());
        }
    }

    if let Some(hdr) = unmapped.header() {
        header = header.set_header(hdr.clone());
    }

    Ok(header.build())
}

/// Copies tags from source record to destination record
///
/// Copies all tags from `src` to `dest`, with special handling:
/// - PG tag is only copied if not already present in dest
/// - Tags in `tag_info.remove` are skipped
///
/// # Arguments
///
/// * `src` - Source record to copy tags from
/// * `dest` - Destination record to copy tags to
/// * `tag_info` - Information about which tags to skip
///
/// # Returns
///
/// `Ok(())` on success
///
/// # Errors
///
/// Currently never returns an error, but signature allows for future error handling
pub fn copy_tags(src: &RecordBuf, dest: &mut Data, tag_info: &TagInfo) -> Result<()> {
    for (tag, value) in src.data().iter() {
        let tag_str = format!("{}{}", tag.as_ref()[0] as char, tag.as_ref()[1] as char);

        if tag_str == "PG" && dest.get(&tag).is_some() {
            continue;
        }

        if tag_info.remove.contains(&tag_str) {
            continue;
        }

        dest.insert(tag, value.clone());
    }
    Ok(())
}

/// Merges tags from unmapped template into mapped template
///
/// This is the core merge function that:
/// 1. Removes specified tags from mapped reads
/// 2. Copies tags from unmapped to mapped reads
/// 3. Applies reverse/revcomp transformations for negative strand
/// 4. Transfers QC pass/fail flags
///
/// Handles both paired-end and single-end data, matching R1 to R1 and R2 to R2.
/// Secondary and supplementary alignments are also processed.
///
/// # Arguments
///
/// * `unmapped` - Template containing unmapped reads with original tags
/// * `mapped` - Template containing mapped reads (will be modified in place)
/// * `tag_info` - Information about which tags to manipulate
///
/// # Returns
///
/// `Ok(())` on success
///
/// # Errors
///
/// Returns an error if tag copying fails (currently never happens)
pub fn merge(unmapped: &Template, mapped: &mut Template, tag_info: &TagInfo) -> Result<()> {
    // Fix mate info first
    mapped.fix_mate_info()?;

    // First, remove tags from mapped reads
    for i in 0..mapped.read_count() {
        for tag_str in &tag_info.remove {
            if tag_str.len() == 2 {
                let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
                let _ = mapped.records[i].data_mut().remove(&tag);
            }
        }
    }

    // Process each primary unmapped read
    for u in unmapped.primary_reads() {
        let is_unpaired = !u.flags().is_segmented();
        let is_first = u.flags().is_first_segment();

        // Find corresponding mapped reads (R1 or R2)
        let mapped_indices: Vec<usize> = if is_unpaired || is_first {
            mapped
                .records
                .iter()
                .enumerate()
                .filter(|(_, r)| !r.flags().is_segmented() || r.flags().is_first_segment())
                .map(|(i, _)| i)
                .collect()
        } else {
            mapped
                .records
                .iter()
                .enumerate()
                .filter(|(_, r)| r.flags().is_segmented() && !r.flags().is_first_segment())
                .map(|(i, _)| i)
                .collect()
        };

        let has_pos_reads =
            mapped_indices.iter().any(|&i| !mapped.records[i].flags().is_reverse_complemented());
        let has_neg_reads =
            mapped_indices.iter().any(|&i| mapped.records[i].flags().is_reverse_complemented());

        // Copy tags to positive strand reads
        if has_pos_reads {
            for &i in &mapped_indices {
                if !mapped.records[i].flags().is_reverse_complemented() {
                    copy_tags(u, mapped.records[i].data_mut(), tag_info)?;
                }
            }
        }

        // Copy tags to negative strand reads (with reverse/revcomp)
        if has_neg_reads {
            let mut u_for_neg = u.clone();

            if tag_info.has_revs_or_revcomps() {
                for (tag, value) in u.data().iter() {
                    let tag_str = format!("{}{}", tag.as_ref()[0] as char, tag.as_ref()[1] as char);

                    if tag_info.reverse.contains(&tag_str) {
                        let new_val = reverse_buf_value(value);
                        u_for_neg.data_mut().insert(tag, new_val);
                    } else if tag_info.revcomp.contains(&tag_str) {
                        let new_val = revcomp_buf_value(value);
                        u_for_neg.data_mut().insert(tag, new_val);
                    }
                }
            }

            for &i in &mapped_indices {
                if mapped.records[i].flags().is_reverse_complemented() {
                    copy_tags(&u_for_neg, mapped.records[i].data_mut(), tag_info)?;
                }
            }
        }

        // Transfer QC pass/fail flag
        let is_qc_fail = u.flags().is_qc_fail();
        for &i in &mapped_indices {
            let mut flags = mapped.records[i].flags();
            if is_qc_fail {
                flags.insert(noodles::sam::alignment::record::Flags::QC_FAIL);
            } else {
                flags.remove(noodles::sam::alignment::record::Flags::QC_FAIL);
            }
            // Apply the modified flags back to the record
            *mapped.records[i].flags_mut() = flags;
        }
    }

    // Normalize alignment integer tags (AS, XS) to Int8 to match fgbio's encoding.
    // When reading from SAM format, noodles may parse these tags with different types
    // (e.g., Int32 or UInt8). fgbio consistently uses signed int8 (SAM type code 'c').
    let as_tag = Tag::new(b'A', b'S');
    let xs_tag = Tag::new(b'X', b'S');
    for i in 0..mapped.read_count() {
        // Normalize AS tag
        if let Some(value) = mapped.records[i].data().get(&as_tag) {
            if let Some(int8_value) = buf_value_to_smallest_signed_int(value) {
                mapped.records[i].data_mut().insert(as_tag, int8_value);
            }
        }
        // Normalize XS tag
        if let Some(value) = mapped.records[i].data().get(&xs_tag) {
            if let Some(int8_value) = buf_value_to_smallest_signed_int(value) {
                mapped.records[i].data_mut().insert(xs_tag, int8_value);
            }
        }
    }

    Ok(())
}

impl Zipper {
    /// Process templates in single-threaded mode
    fn process_singlethreaded<W, U, M>(
        &self,
        unmapped_iter: U,
        mut mapped_iter: M,
        output_header: &Header,
        writer: &mut W,
        tag_info: &TagInfo,
    ) -> Result<u64>
    where
        W: AlignmentWrite,
        U: Iterator<Item = Result<Template>>,
        M: Iterator<Item = Result<Template>>,
    {
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
        let mut mapped_peek: Option<Template> = None;
        let mut templates_not_in_mapped_bam: u64 = 0;

        for unmapped_result in unmapped_iter {
            let unmapped_template = unmapped_result?;

            if mapped_peek.is_none() {
                mapped_peek = mapped_iter.next().transpose()?;
            }

            if let Some(ref mut mapped_template) = mapped_peek {
                if mapped_template.name == unmapped_template.name {
                    merge(&unmapped_template, mapped_template, tag_info)?;
                    for rec in mapped_template.all_reads() {
                        writer.write_alignment_record(output_header, rec)?;
                        progress.log_if_needed(1);
                    }
                    mapped_peek = None;
                } else {
                    // Unmapped read with no corresponding mapped read (orphan unmapped)
                    debug!(
                        "Found unmapped read with no corresponding mapped read: {}",
                        String::from_utf8_lossy(&unmapped_template.name)
                    );
                    if self.exclude_missing_reads {
                        templates_not_in_mapped_bam += 1;
                    } else {
                        for rec in unmapped_template.all_reads() {
                            writer.write_alignment_record(output_header, rec)?;
                            progress.log_if_needed(1);
                        }
                    }
                }
            } else {
                // Unmapped read with no more mapped reads (orphan unmapped)
                debug!(
                    "Found unmapped read with no corresponding mapped read: {}",
                    String::from_utf8_lossy(&unmapped_template.name)
                );
                if self.exclude_missing_reads {
                    templates_not_in_mapped_bam += 1;
                } else {
                    for rec in unmapped_template.all_reads() {
                        writer.write_alignment_record(output_header, rec)?;
                        progress.log_if_needed(1);
                    }
                }
            }
        }

        progress.log_final();

        // Check for leftover mapped reads - this is an error
        if let Some(remaining) = mapped_peek {
            anyhow::bail!(
                "Error: processed all unmapped reads but there are mapped reads remaining. \
                 Found template '{}'. Please ensure the unmapped and mapped reads have the same set \
                 of read names in the same order, and reads with the same name are consecutive (grouped) \
                 in each input.",
                String::from_utf8_lossy(&remaining.name)
            );
        }
        if let Some(remaining) = mapped_iter.next() {
            let remaining = remaining?;
            anyhow::bail!(
                "Error: processed all unmapped reads but there are mapped reads remaining. \
                 Found template '{}'. Please ensure the unmapped and mapped reads have the same set \
                 of read names in the same order, and reads with the same name are consecutive (grouped) \
                 in each input.",
                String::from_utf8_lossy(&remaining.name)
            );
        }

        if self.exclude_missing_reads && templates_not_in_mapped_bam > 0 {
            info!(
                "Excluded {templates_not_in_mapped_bam} templates that were not present in the aligned BAM."
            );
        }

        Ok(progress.count())
    }

    /// Process templates with multi-threaded BGZF compression.
    ///
    /// This method uses noodles' multi-threaded writer to parallelize compression
    /// (the bottleneck) while processing templates sequentially.
    fn process_with_parallel_writer<U, M>(
        &self,
        unmapped_iter: U,
        mut mapped_iter: M,
        output_header: &Header,
        tag_info: &TagInfo,
    ) -> Result<u64>
    where
        U: Iterator<Item = Result<Template>>,
        M: Iterator<Item = Result<Template>>,
    {
        // Create multi-threaded BAM writer
        let mut writer = create_bam_writer(
            &self.output,
            output_header,
            self.threads,
            self.compression.compression_level,
        )?;

        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);
        let mut mapped_peek: Option<Template> = None;
        let mut templates_not_in_mapped_bam: u64 = 0;

        for unmapped_result in unmapped_iter {
            let unmapped_template = unmapped_result?;

            // Advance mapped iterator if needed
            if mapped_peek.is_none() {
                mapped_peek = mapped_iter.next().transpose()?;
            }

            // Check for matching mapped template
            if let Some(ref mut mapped_template) = mapped_peek {
                if mapped_template.name == unmapped_template.name {
                    // Merge tags from unmapped to mapped (no cloning needed!)
                    merge(&unmapped_template, mapped_template, tag_info)?;

                    // Write merged records
                    for rec in mapped_template.all_reads() {
                        writer.write_alignment_record(output_header, rec)?;
                        progress.log_if_needed(1);
                    }
                    mapped_peek = None;
                } else {
                    // Unmapped read with no corresponding mapped read (orphan unmapped)
                    debug!(
                        "Found unmapped read with no corresponding mapped read: {}",
                        String::from_utf8_lossy(&unmapped_template.name)
                    );
                    if self.exclude_missing_reads {
                        templates_not_in_mapped_bam += 1;
                    } else {
                        for rec in unmapped_template.all_reads() {
                            writer.write_alignment_record(output_header, rec)?;
                            progress.log_if_needed(1);
                        }
                    }
                }
            } else {
                // No more mapped reads - write unmapped reads
                debug!(
                    "Found unmapped read with no corresponding mapped read: {}",
                    String::from_utf8_lossy(&unmapped_template.name)
                );
                if self.exclude_missing_reads {
                    templates_not_in_mapped_bam += 1;
                } else {
                    for rec in unmapped_template.all_reads() {
                        writer.write_alignment_record(output_header, rec)?;
                        progress.log_if_needed(1);
                    }
                }
            }
        }

        progress.log_final();

        // Check for leftover mapped reads - this is an error
        if let Some(remaining) = mapped_peek {
            // Finish writer before returning error
            writer.into_inner().finish()?;
            anyhow::bail!(
                "Error: processed all unmapped reads but there are mapped reads remaining. \
                 Found template '{}'. Please ensure the unmapped and mapped reads have the same set \
                 of read names in the same order, and reads with the same name are consecutive (grouped) \
                 in each input.",
                String::from_utf8_lossy(&remaining.name)
            );
        }
        if let Some(remaining) = mapped_iter.next() {
            let remaining = remaining?;
            // Finish writer before returning error
            writer.into_inner().finish()?;
            anyhow::bail!(
                "Error: processed all unmapped reads but there are mapped reads remaining. \
                 Found template '{}'. Please ensure the unmapped and mapped reads have the same set \
                 of read names in the same order, and reads with the same name are consecutive (grouped) \
                 in each input.",
                String::from_utf8_lossy(&remaining.name)
            );
        }

        if self.exclude_missing_reads && templates_not_in_mapped_bam > 0 {
            info!(
                "Excluded {templates_not_in_mapped_bam} templates that were not present in the aligned BAM."
            );
        }

        // Finish writing and flush all pending blocks
        writer.into_inner().finish()?;

        Ok(progress.count())
    }
}

impl Command for Zipper {
    /// Executes the `zipper` command
    ///
    /// Main workflow:
    /// 1. Opens and validates input BAM files
    /// 2. Checks queryname sorting
    /// 3. Builds output header from both inputs and reference dictionary
    /// 4. Iterates through templates from both files in parallel
    /// 5. Merges tags from unmapped to mapped reads
    /// 6. Writes merged records to output
    ///
    /// # Returns
    ///
    /// `Ok(())` on successful completion
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Input files cannot be opened
    /// - Files are not properly queryname sorted/grouped
    /// - Reference dictionary is missing or invalid
    /// - Input files have different sets of read names
    /// - I/O errors occur during reading or writing
    fn execute(&self, command_line: &str) -> Result<()> {
        info!("Starting zipper");

        let timer = OperationTimer::new("Zipping BAMs");

        // Validate input files exist
        validate_file_exists(&self.unmapped, "unmapped BAM file")?;
        validate_file_exists(&self.reference, "reference FASTA file")?;
        let dict_path = self.reference.with_extension("dict");
        anyhow::ensure!(
            dict_path.exists(),
            "Reference dictionary file does not exist: {}",
            dict_path.display()
        );

        let (mut unmapped_reader, unmapped_header) = create_bam_reader(&self.unmapped, 1)?;

        // For stdin, use BatchedSamReader that starts at 64MB and grows based on observed usage.
        // This handles bwa mem's bursty output pattern (450-750MB SAM text per batch).
        // The -K parameter (bwa_chunk_size) is used for batch tracking.
        let mapped_reader: Box<dyn BufRead + Send> = if self.input.to_str() == Some("-") {
            info!("Reading SAM from stdin with adaptive buffer (bwa -K {})", self.bwa_chunk_size);
            Box::new(BatchedSamReader::new(std::io::stdin(), self.bwa_chunk_size))
        } else {
            Box::new(BufReader::with_capacity(
                256 * 1024,
                std::fs::File::open(&self.input).context("Failed to open mapped SAM")?,
            ))
        };
        let mut mapped_sam_reader = noodles::sam::io::Reader::new(mapped_reader);
        let mapped_header = mapped_sam_reader.read_header()?;

        check_sort(&unmapped_header, &self.unmapped, "unmapped");
        check_sort(&mapped_header, &self.input, "mapped");

        let output_header = build_output_header(&unmapped_header, &mapped_header, &self.reference)?;

        // Add @PG record with PP chaining
        let output_header = fgumi_lib::header::add_pg_record(
            output_header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        let tag_info = TagInfo::new(
            self.tags_to_remove.clone(),
            self.tags_to_reverse.clone(),
            self.tags_to_revcomp.clone(),
        );

        if !tag_info.remove.is_empty() {
            info!("Tags for removal: {:?}", tag_info.remove);
        }
        if !tag_info.reverse.is_empty() {
            info!("Tags being reversed: {:?}", tag_info.reverse);
        }
        if !tag_info.revcomp.is_empty() {
            info!("Tags being reverse complemented: {:?}", tag_info.revcomp);
        }

        // Create async unmapped reader - spawn thread to read ahead
        let (unmapped_tx, unmapped_rx) =
            std::sync::mpsc::sync_channel::<Result<Template>>(self.buffer);
        let unmapped_header_for_reader = unmapped_header.clone();
        std::thread::spawn(move || {
            let unmapped_iter = TemplateIterator::new(
                unmapped_reader
                    .record_bufs(&unmapped_header_for_reader)
                    .map(|r| r.map_err(anyhow::Error::from)),
            );
            for template in unmapped_iter {
                if unmapped_tx.send(template).is_err() {
                    break; // Receiver dropped, main thread done
                }
            }
        });
        let unmapped_iter = std::iter::from_fn(move || unmapped_rx.recv().ok());

        // Create async mapped reader - spawn thread to read ahead and prevent backpressure on stdin
        // (e.g., when reading piped output from bwa)
        let (mapped_tx, mapped_rx) = std::sync::mpsc::sync_channel::<Result<Template>>(self.buffer);
        let mapped_header_for_reader = mapped_header.clone();
        std::thread::spawn(move || {
            let mapped_iter = TemplateIterator::new(
                mapped_sam_reader
                    .record_bufs(&mapped_header_for_reader)
                    .map(|r| r.map_err(anyhow::Error::from)),
            );
            for template in mapped_iter {
                if mapped_tx.send(template).is_err() {
                    break; // Receiver dropped, main thread done
                }
            }
        });
        let mapped_iter = std::iter::from_fn(move || mapped_rx.recv().ok());

        let total_records = if self.threads > 1 {
            // Multi-threaded: sequential merge with parallel BGZF compression
            // This is more efficient than parallel merge with template cloning
            self.process_with_parallel_writer(
                unmapped_iter,
                mapped_iter,
                &output_header,
                &tag_info,
            )?
        } else {
            // Single-threaded: simple sequential processing
            let output_writer: Box<dyn Write> = if self.output.to_str() == Some("-") {
                Box::new(std::io::stdout().lock())
            } else {
                Box::new(
                    std::fs::File::create(&self.output).context("Failed to create output BAM")?,
                )
            };

            let mut writer = noodles::bam::io::Writer::new(output_writer);
            writer.write_header(&output_header)?;

            self.process_singlethreaded(
                unmapped_iter,
                mapped_iter,
                &output_header,
                &mut writer,
                &tag_info,
            )?
        };

        info!("zipper completed successfully");
        timer.log_completion(total_records);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use bstr::ByteSlice;
    use fgumi_lib::sam::builder::{
        MAPPED_PG_ID, REFERENCE_LENGTH, RecordBuilder, SamBuilder, create_ref_dict,
    };
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
    use std::collections::HashMap;
    use tempfile::TempDir;

    /// Runs `zipper-bams` and returns the output records
    ///
    /// Helper function that:
    /// 1. Writes unmapped and mapped builders to temporary BAM files
    /// 2. Creates a temporary reference dictionary
    /// 3. Executes `zipper-bams`
    /// 4. Reads and returns all output records
    ///
    /// # Arguments
    ///
    /// * `unmapped` - Builder containing unmapped reads
    /// * `mapped` - Builder containing mapped reads
    /// * `tags_to_remove` - Tags to remove from mapped reads
    /// * `tags_to_reverse` - Tags to reverse for negative strand
    /// * `tags_to_revcomp` - Tags to revcomp for negative strand
    ///
    /// # Returns
    ///
    /// Vector of all records written to output BAM
    ///
    /// # Errors
    ///
    /// Returns an error if any step fails
    fn run_zipper(
        unmapped: &SamBuilder,
        mapped: &SamBuilder,
        tags_to_remove: Vec<String>,
        tags_to_reverse: Vec<String>,
        tags_to_revcomp: Vec<String>,
    ) -> Result<Vec<RecordBuf>> {
        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");
        let dict_path = create_ref_dict(&dir, "chr1", REFERENCE_LENGTH)?;

        unmapped.write(&unmapped_path)?;
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path.clone(),
            unmapped: unmapped_path.clone(),
            reference: dict_path,
            output: output_path.clone(),
            tags_to_remove,
            tags_to_reverse,
            tags_to_revcomp,
            buffer: 5000,
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads: false,
        };

        zipper.execute("test")?;

        // Read output
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&output_path)?;
        let header = reader.read_header()?;

        let mut records = Vec::new();
        let mut record = RecordBuf::default();
        loop {
            match reader.read_record_buf(&header, &mut record) {
                Ok(0) => break,
                Ok(_) => {
                    records.push(record.clone());
                }
                Err(e) => return Err(e.into()),
            }
        }

        Ok(records)
    }

    /// Tests basic tag merging from unmapped to mapped reads
    ///
    /// Verifies that:
    /// - Tags from unmapped reads are copied to mapped reads
    /// - Tags from mapped reads are preserved
    /// - All records are present in output
    /// - MC and MQ tags are properly set for paired reads
    #[test]
    fn test_basic_merge() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped pairs with tags
        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));
        attrs1.insert("xy", BufValue::from(1234i32));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs1);

        let mut attrs2 = HashMap::new();
        attrs2.insert("RX", BufValue::from("GGTA".to_string()));
        attrs2.insert("xy", BufValue::from(4567i32));
        unmapped.add_pair_with_attrs("q2", None, None, true, true, &attrs2);

        // Add mapped pairs with different tags
        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs1);

        let mut mapped_attrs2 = HashMap::new();
        mapped_attrs2.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs2.insert("AS", BufValue::from(16i32));
        mapped.add_pair_with_attrs("q2", Some(500), Some(700), true, true, &mapped_attrs2);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 4);

        // Check that tags from unmapped BAM were copied
        for rec in &records {
            // Should have RX tag from unmapped
            assert!(rec.data().get(&Tag::new(b'R', b'X')).is_some());

            // Should have xy tag from unmapped
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());

            // Should have AS tag from mapped
            assert!(rec.data().get(&Tag::new(b'A', b'S')).is_some());

            // Should have PG tag from mapped
            assert!(rec.data().get(&Tag::new(b'P', b'G')).is_some());

            // Should have MC (mate CIGAR) and MQ (mate mapping quality) tags
            assert!(rec.data().get(&Tag::new(b'M', b'C')).is_some(), "MC tag should be present");
            assert!(rec.data().get(&Tag::new(b'M', b'Q')).is_some(), "MQ tag should be present");
        }

        Ok(())
    }

    /// Tests handling of unpaired (single-end) reads
    ///
    /// Verifies that single-end reads are properly merged without
    /// attempting to match R1/R2.
    #[test]
    fn test_unpaired_reads() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));
        attrs1.insert("xy", BufValue::from(1234i32));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs1);

        let mut attrs2 = HashMap::new();
        attrs2.insert("RX", BufValue::from("GGTA".to_string()));
        attrs2.insert("xy", BufValue::from(4567i32));
        unmapped.add_frag_with_attrs("q2", None, true, &attrs2);

        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs1);

        let mut mapped_attrs2 = HashMap::new();
        mapped_attrs2.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs2.insert("AS", BufValue::from(16i32));
        mapped.add_frag_with_attrs("q2", Some(500), false, &mapped_attrs2);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        for rec in &records {
            // Should have tags from both unmapped and mapped
            assert!(rec.data().get(&Tag::new(b'R', b'X')).is_some());
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());
            assert!(rec.data().get(&Tag::new(b'A', b'S')).is_some());
        }

        Ok(())
    }

    /// Tests tag removal, reversal, and reverse complementation
    ///
    /// Verifies that:
    /// - Specified tags are removed from mapped reads
    /// - Array and string tags are reversed for negative strand
    /// - Sequence tags are reverse complemented for negative strand
    /// - Positive strand reads are unchanged
    #[test]
    fn test_tag_removal_reverse_revcomp() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("n1", BufValue::from(vec![1i16, 2, 3, 4, 5]));
        attrs.insert("n2", BufValue::from(vec![2i16, 3, 4, 5, 6]));
        attrs.insert("s1", BufValue::from("abcde".to_string()));
        attrs.insert("s2", BufValue::from("vwxyz".to_string()));
        attrs.insert("s3", BufValue::from("AGAGG".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, false, &mapped_attrs);

        let records = run_zipper(
            &unmapped,
            &mapped,
            vec!["AS".to_string()],
            vec!["n1".to_string(), "n2".to_string(), "s1".to_string()],
            vec!["s3".to_string()],
        )?;

        assert_eq!(records.len(), 2);

        for rec in &records {
            // AS tag should be removed
            assert!(rec.data().get(&Tag::new(b'A', b'S')).is_none());

            if rec.flags().is_first_segment() {
                // R1 is positive strand - no reversing/revcomping
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(*vals, vec![1i16, 2, 3, 4, 5]);
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(s.to_string(), "abcde");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(s.to_string(), "AGAGG");
                }
            } else {
                // R2 is negative strand - should be reversed/revcomped
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(*vals, vec![5i16, 4, 3, 2, 1]);
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(s.to_string(), "edcba");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(s.to_string(), "CCTCT");
                }
            }

            // s2 should not be changed
            if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'2')) {
                assert_eq!(s.to_string(), "vwxyz");
            }
        }

        Ok(())
    }

    /// Tests the "Consensus" named tag set expansion
    ///
    /// Verifies that specifying "Consensus" in `tags_to_reverse` or `tags_to_revcomp`
    /// properly expands to all per-base consensus tags and applies the correct
    /// transformations.
    #[test]
    fn test_consensus_tag_set() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("aD", BufValue::from("AAAGG".to_string()));
        attrs.insert("bD", BufValue::from("AAAGC".to_string()));
        attrs.insert("ad", BufValue::from(vec![3i16, 3, 4, 4, 2]));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs1.clone());

        // Add supplementary record
        let mut supp_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SUPPLEMENTARY | Flags::REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        for (tag_str, value) in &mapped_attrs1 {
            let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
            supp_rec.data_mut().insert(tag, value.clone());
        }
        mapped.push_record(supp_rec);

        let records = run_zipper(
            &unmapped,
            &mapped,
            vec![],
            vec!["Consensus".to_string()],
            vec!["Consensus".to_string()],
        )?;

        assert_eq!(records.len(), 2);

        for rec in &records {
            if rec.flags().is_reverse_complemented() {
                // Negative strand - should be reversed/revcomped
                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'a', b'D')) {
                    assert_eq!(s.to_string(), "CCTTT");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'b', b'D')) {
                    assert_eq!(s.to_string(), "GCTTT");
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'a', b'd'))
                {
                    assert_eq!(*vals, vec![2i16, 4, 4, 3, 3]);
                }
            } else {
                // Positive strand - no changes
                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'a', b'D')) {
                    assert_eq!(s.to_string(), "AAAGG");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b'b', b'D')) {
                    assert_eq!(s.to_string(), "AAAGC");
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'a', b'd'))
                {
                    assert_eq!(*vals, vec![3i16, 3, 4, 4, 2]);
                }
            }
        }

        Ok(())
    }

    /// Tests handling of reads that are only in the unmapped BAM
    ///
    /// Verifies that if a read has no corresponding mapped read, it is
    /// written to the output unchanged (still unmapped).
    #[test]
    fn test_unmapped_only_reads() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));
        attrs1.insert("xy", BufValue::from(1234i32));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs1);

        let mut attrs2 = HashMap::new();
        attrs2.insert("RX", BufValue::from("GATA".to_string()));
        attrs2.insert("xy", BufValue::from(3456i32));
        unmapped.add_frag_with_attrs("q2", None, true, &attrs2);

        let mut attrs3 = HashMap::new();
        attrs3.insert("RX", BufValue::from("GGCG".to_string()));
        attrs3.insert("xy", BufValue::from(5678i32));
        unmapped.add_frag_with_attrs("q3", None, true, &attrs3);

        let mut mapped_attrs1 = HashMap::new();
        mapped_attrs1.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs1.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs1);

        let mut mapped_attrs3 = HashMap::new();
        mapped_attrs3.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs3.insert("AS", BufValue::from(77i32));
        mapped.add_frag_with_attrs("q3", Some(200), false, &mapped_attrs3);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 3);

        let names: Vec<String> = records
            .iter()
            .map(|r| String::from_utf8(r.name().unwrap().as_bytes().to_vec()).unwrap())
            .collect();

        assert_eq!(names, vec!["q1", "q2", "q3"]);

        // q2 should be unmapped
        let q2 = records
            .iter()
            .find(|r| String::from_utf8(r.name().unwrap().as_bytes().to_vec()).unwrap() == "q2")
            .unwrap();
        assert!(q2.flags().is_unmapped());

        Ok(())
    }

    /// Tests `TagInfo` expansion of named tag sets
    ///
    /// Verifies that the "Consensus" tag set is properly expanded into
    /// its constituent per-base tags.
    #[test]
    fn test_tag_info_expansion() {
        let tag_info =
            TagInfo::new(vec![], vec!["Consensus".to_string()], vec!["Consensus".to_string()]);

        // Check that "Consensus" was expanded
        assert!(tag_info.reverse.contains("ad"));
        assert!(tag_info.reverse.contains("ae"));
        assert!(tag_info.reverse.contains("bd"));
        assert!(tag_info.reverse.contains("be"));
        assert!(tag_info.reverse.contains("cd"));

        assert!(tag_info.revcomp.contains("aD"));
        assert!(tag_info.revcomp.contains("bD"));
        assert!(tag_info.revcomp.contains("cD"));
    }

    #[test]
    fn test_mate_info_fixing() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        // Find R1 and R2
        let r1 = records.iter().find(|r| r.flags().is_first_segment()).unwrap();
        let r2 = records.iter().find(|r| !r.flags().is_first_segment()).unwrap();

        // Check that MQ tags are set
        let r1_mq = r1.data().get(&Tag::new(b'M', b'Q'));
        let r2_mq = r2.data().get(&Tag::new(b'M', b'Q'));
        assert!(r1_mq.is_some(), "R1 should have MQ tag");
        assert!(r2_mq.is_some(), "R2 should have MQ tag");

        // Check that MC tags are set
        let r1_mc = r1.data().get(&Tag::new(b'M', b'C'));
        let r2_mc = r2.data().get(&Tag::new(b'M', b'C'));
        assert!(r1_mc.is_some(), "R1 should have MC tag");
        assert!(r2_mc.is_some(), "R2 should have MC tag");

        Ok(())
    }

    /// Tests QC flag transfer from unmapped to mapped reads
    ///
    /// Verifies that:
    /// - `QC_FAIL` flag is transferred from unmapped reads to mapped reads
    /// - `QC_PASS` (absence of flag) is also properly transferred
    /// - Both R1 and R2 get the correct flags
    #[test]
    fn test_qc_flag_transfer() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped pair where R1 passes QC but R2 fails QC
        let mut attrs1 = HashMap::new();
        attrs1.insert("RX", BufValue::from("ACGT".to_string()));

        // Create R1 (passes QC - no QC_FAIL flag)
        let r1_unmapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::UNMAPPED)
            .tag("RX", "ACGT")
            .build();

        // Create R2 (fails QC - has QC_FAIL flag)
        let r2_unmapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::UNMAPPED | Flags::QC_FAIL)
            .tag("RX", "ACGT")
            .build();

        unmapped.push_record(r1_unmapped);
        unmapped.push_record(r2_unmapped);

        // Add corresponding mapped pair (neither has QC_FAIL initially)
        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        // Find R1 and R2 in output
        let r1 = records.iter().find(|r| r.flags().is_first_segment()).unwrap();
        let r2 = records.iter().find(|r| !r.flags().is_first_segment()).unwrap();

        // R1 should NOT have QC_FAIL flag (it passed QC in unmapped)
        assert!(
            !r1.flags().is_qc_fail(),
            "R1 should not have QC_FAIL flag (passed QC in unmapped)"
        );

        // R2 SHOULD have QC_FAIL flag (it failed QC in unmapped)
        assert!(r2.flags().is_qc_fail(), "R2 should have QC_FAIL flag (failed QC in unmapped)");

        Ok(())
    }

    /// Tests QC flag removal when unmapped read passes QC but mapped read had failed
    ///
    /// Verifies that if a mapped read has `QC_FAIL` but the unmapped read doesn't,
    /// the flag is removed from the mapped read.
    #[test]
    fn test_qc_flag_removal() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped read that passes QC
        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        // Add mapped read with QC_FAIL flag set
        let mapped_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::QC_FAIL)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        mapped.push_record(mapped_rec);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 1);

        // The QC_FAIL flag should have been removed
        assert!(
            !records[0].flags().is_qc_fail(),
            "QC_FAIL flag should be removed when unmapped read passes QC"
        );

        Ok(())
    }

    /// Tests tag removal, reversal, and reverse complementation on secondary and supplementary alignments
    ///
    /// Verifies that:
    /// - Tag operations work correctly on both secondary and supplementary records
    /// - Multiple alignments (primary, secondary, supplementary) all get processed
    /// - Negative strand transformations apply to secondary/supplementary alignments
    #[test]
    #[allow(clippy::too_many_lines)]
    fn test_tag_operations_on_secondary_and_supplementary() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("n1", BufValue::from(vec![1i16, 2, 3, 4, 5]));
        attrs.insert("n2", BufValue::from(vec![2i16, 3, 4, 5, 6]));
        attrs.insert("s1", BufValue::from("abcde".to_string()));
        attrs.insert("s2", BufValue::from("vwxyz".to_string()));
        attrs.insert("s3", BufValue::from("AGAGG".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped_attrs.insert("AS", BufValue::from(77i32));

        // Add primary alignment (positive strand)
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs.clone());

        // Add secondary alignment (negative strand)
        let mut secondary_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SECONDARY | Flags::REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        for (tag_str, value) in &mapped_attrs {
            let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
            secondary_rec.data_mut().insert(tag, value.clone());
        }
        mapped.push_record(secondary_rec);

        // Add supplementary alignment (negative strand)
        let mut supp_rec = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SUPPLEMENTARY | Flags::REVERSE_COMPLEMENTED)
            .reference_sequence_id(0)
            .alignment_start(300)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .build();
        for (tag_str, value) in &mapped_attrs {
            let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
            supp_rec.data_mut().insert(tag, value.clone());
        }
        mapped.push_record(supp_rec);

        let records = run_zipper(
            &unmapped,
            &mapped,
            vec!["AS".to_string()],
            vec!["n1".to_string(), "n2".to_string(), "s1".to_string()],
            vec!["s3".to_string()],
        )?;

        assert_eq!(records.len(), 3, "Should have 3 records (primary, secondary, supplementary)");

        let secondary_and_supp_count = records
            .iter()
            .filter(|r| r.flags().is_secondary() || r.flags().is_supplementary())
            .count();
        assert_eq!(secondary_and_supp_count, 2, "Should have 2 secondary/supplementary records");

        let neg_strand_count =
            records.iter().filter(|r| r.flags().is_reverse_complemented()).count();
        assert_eq!(neg_strand_count, 2, "Should have 2 negative strand records");

        // Check each record
        for rec in &records {
            // AS tag should be removed from all records
            assert!(rec.data().get(&Tag::new(b'A', b'S')).is_none(), "AS tag should be removed");

            if rec.flags().is_reverse_complemented() {
                // Negative strand - should be reversed/revcomped
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(
                        *vals,
                        vec![5i16, 4, 3, 2, 1],
                        "n1 should be reversed for negative strand"
                    );
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'2'))
                {
                    assert_eq!(
                        *vals,
                        vec![6i16, 5, 4, 3, 2],
                        "n2 should be reversed for negative strand"
                    );
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(s.to_string(), "edcba", "s1 should be reversed for negative strand");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'2')) {
                    assert_eq!(s.to_string(), "vwxyz", "s2 should not be changed");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(
                        s.to_string(),
                        "CCTCT",
                        "s3 should be revcomped for negative strand"
                    );
                }
            } else {
                // Positive strand - no changes
                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'1'))
                {
                    assert_eq!(
                        *vals,
                        vec![1i16, 2, 3, 4, 5],
                        "n1 should be unchanged for positive strand"
                    );
                }

                if let Some(BufValue::Array(
                    noodles::sam::alignment::record_buf::data::field::value::Array::Int16(vals),
                )) = rec.data().get(&Tag::new(b'n', b'2'))
                {
                    assert_eq!(
                        *vals,
                        vec![2i16, 3, 4, 5, 6],
                        "n2 should be unchanged for positive strand"
                    );
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'1')) {
                    assert_eq!(
                        s.to_string(),
                        "abcde",
                        "s1 should be unchanged for positive strand"
                    );
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'2')) {
                    assert_eq!(s.to_string(), "vwxyz", "s2 should not be changed");
                }

                if let Some(BufValue::String(s)) = rec.data().get(&Tag::new(b's', b'3')) {
                    assert_eq!(
                        s.to_string(),
                        "AGAGG",
                        "s3 should be unchanged for positive strand"
                    );
                }
            }
        }

        Ok(())
    }

    /// Tests multi-threaded processing
    ///
    /// Verifies that multi-threaded mode produces the same results as single-threaded
    #[test]
    fn test_multithreaded_processing() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add multiple templates to process
        for i in 0..10 {
            let name = format!("q{i}");
            let mut attrs = HashMap::new();
            attrs.insert("RX", BufValue::from(format!("ACG{i}")));
            attrs.insert("xy", BufValue::from(1000 + i as i32));
            unmapped.add_pair_with_attrs(&name, None, None, true, true, &attrs);

            let mut mapped_attrs = HashMap::new();
            mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
            mapped_attrs.insert("AS", BufValue::from(77i32));
            mapped.add_pair_with_attrs(
                &name,
                Some(100 + i * 100),
                Some(200 + i * 100),
                true,
                true,
                &mapped_attrs,
            );
        }

        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");
        let dict_path = create_ref_dict(&dir, "chr1", REFERENCE_LENGTH)?;

        unmapped.write(&unmapped_path)?;
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path.clone(),
            unmapped: unmapped_path.clone(),
            reference: dict_path,
            output: output_path.clone(),
            tags_to_remove: vec![],
            tags_to_reverse: vec![],
            tags_to_revcomp: vec![],
            buffer: 5000,
            threads: 4,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads: false,
        };

        zipper.execute("test")?;

        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&output_path)?;
        let header = reader.read_header()?;

        let mut records = Vec::new();
        let mut record = RecordBuf::default();
        loop {
            match reader.read_record_buf(&header, &mut record) {
                Ok(0) => break,
                Ok(_) => {
                    records.push(record.clone());
                }
                Err(e) => return Err(e.into()),
            }
        }

        assert_eq!(records.len(), 20);

        // Verify all records have both unmapped and mapped tags
        for rec in &records {
            assert!(rec.data().get(&Tag::new(b'R', b'X')).is_some());
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());
            assert!(rec.data().get(&Tag::new(b'A', b'S')).is_some());
        }

        Ok(())
    }

    /// Tests empty unmapped input with non-empty mapped
    ///
    /// Verifies that having mapped reads remaining after unmapped is exhausted is an error
    /// (matching fgbio's `ZipperBams` validation)
    #[test]
    fn test_empty_unmapped() -> Result<()> {
        let unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs);

        let result = run_zipper(&unmapped, &mapped, vec![], vec![], vec![]);

        // With the leftover mapped reads validation (matching fgbio), having mapped
        // reads remaining after unmapped is exhausted should be an error
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("mapped reads remaining"),
            "Expected error about leftover mapped reads, got: {err_msg}"
        );

        Ok(())
    }

    /// Tests `copy_tags` function directly
    ///
    /// Verifies the tag copying logic including PG tag special handling
    #[test]
    fn test_copy_tags_with_pg_present() -> Result<()> {
        let src =
            RecordBuilder::new().sequence("ACGT").tag("PG", "old_pg").tag("XY", 123i32).build();

        let mut dest_data = Data::default();
        dest_data.insert(Tag::new(b'P', b'G'), BufValue::from("new_pg".to_string()));

        let tag_info = TagInfo::new(vec![], vec![], vec![]);
        copy_tags(&src, &mut dest_data, &tag_info)?;

        // PG should not be overwritten
        if let Some(BufValue::String(s)) = dest_data.get(&Tag::new(b'P', b'G')) {
            assert_eq!(s.to_string(), "new_pg");
        }

        // XY should be copied
        assert!(dest_data.get(&Tag::new(b'X', b'Y')).is_some());

        Ok(())
    }

    /// Tests `copy_tags` with tag removal
    ///
    /// Verifies that tags in the remove list are not copied
    #[test]
    fn test_copy_tags_with_removal() -> Result<()> {
        let src = RecordBuilder::new().sequence("ACGT").tag("XY", 123i32).tag("AB", 456i32).build();

        let mut dest_data = Data::default();

        let tag_info = TagInfo::new(vec!["XY".to_string()], vec![], vec![]);
        copy_tags(&src, &mut dest_data, &tag_info)?;

        // XY should not be copied (it's in remove list)
        assert!(dest_data.get(&Tag::new(b'X', b'Y')).is_none());

        // AB should be copied
        assert!(dest_data.get(&Tag::new(b'A', b'B')).is_some());

        Ok(())
    }

    /// Tests merge function with mixed read configurations
    ///
    /// Verifies handling of reads with different strand orientations
    #[test]
    fn test_merge_mixed_strands() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        attrs.insert("xy", BufValue::from(1234i32));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        // Both reads map to positive strand
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 2);

        // All records should have the unmapped tags
        for rec in &records {
            assert!(rec.data().get(&Tag::new(b'R', b'X')).is_some());
            assert!(rec.data().get(&Tag::new(b'x', b'y')).is_some());
        }

        Ok(())
    }

    /// Tests `TagInfo` `has_revs_or_revcomps` method
    ///
    /// Verifies the helper method correctly identifies when reversals/revcomps are configured
    #[test]
    fn test_tag_info_has_revs_or_revcomps() {
        let tag_info1 = TagInfo::new(vec![], vec![], vec![]);
        assert!(!tag_info1.has_revs_or_revcomps());

        let tag_info2 = TagInfo::new(vec![], vec!["n1".to_string()], vec![]);
        assert!(tag_info2.has_revs_or_revcomps());

        let tag_info3 = TagInfo::new(vec![], vec![], vec!["s1".to_string()]);
        assert!(tag_info3.has_revs_or_revcomps());

        let tag_info4 = TagInfo::new(vec![], vec!["n1".to_string()], vec!["s1".to_string()]);
        assert!(tag_info4.has_revs_or_revcomps());
    }

    /// Tests handling of multiple secondary alignments
    ///
    /// Verifies that multiple secondary alignments for the same read are all processed
    #[test]
    fn test_multiple_secondary_alignments() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_frag_with_attrs("q1", None, true, &attrs);

        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_frag_with_attrs("q1", Some(100), true, &mapped_attrs.clone());

        // Add 3 secondary alignments
        for i in 0..3 {
            let mut secondary_rec = RecordBuilder::new()
                .name("q1")
                .flags(Flags::SECONDARY)
                .reference_sequence_id(0)
                .alignment_start(200 + i * 100)
                .cigar("100M")
                .sequence(&"A".repeat(100))
                .qualities(&[30u8; 100])
                .build();
            for (tag_str, value) in &mapped_attrs {
                let tag = Tag::new(tag_str.as_bytes()[0], tag_str.as_bytes()[1]);
                secondary_rec.data_mut().insert(tag, value.clone());
            }
            mapped.push_record(secondary_rec);
        }

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 4, "Should have 1 primary + 3 secondary alignments");

        // All records should have the RX tag
        for rec in &records {
            assert!(rec.data().get(&Tag::new(b'R', b'X')).is_some());
        }

        Ok(())
    }

    /// Tests handling of R2-only mapping (R1 unmapped, R2 mapped)
    ///
    /// Verifies correct handling when only one read of a pair is mapped
    #[test]
    fn test_r2_only_mapping() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);

        // Only R2 is mapped
        let r2_mapped = RecordBuilder::new()
            .name("q1")
            .flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("100M")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .tag("PG", MAPPED_PG_ID)
            .build();
        mapped.push_record(r2_mapped);

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        assert_eq!(records.len(), 1, "Should have only R2 in output");

        let r2 = &records[0];
        assert!(!r2.flags().is_first_segment());
        assert!(r2.data().get(&Tag::new(b'R', b'X')).is_some());

        Ok(())
    }

    /// Helper to run zipper with `exclude_missing_reads` option
    fn run_zipper_with_options(
        unmapped: &SamBuilder,
        mapped: &SamBuilder,
        exclude_missing_reads: bool,
    ) -> Result<Vec<RecordBuf>> {
        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");
        let dict_path = create_ref_dict(&dir, "chr1", REFERENCE_LENGTH)?;

        unmapped.write(&unmapped_path)?;
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path.clone(),
            unmapped: unmapped_path.clone(),
            reference: dict_path,
            output: output_path.clone(),
            tags_to_remove: vec![],
            tags_to_reverse: vec![],
            tags_to_revcomp: vec![],
            buffer: 5000,
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads,
        };

        zipper.execute("test")?;

        // Read output
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&output_path)?;
        let header = reader.read_header()?;

        let mut records = Vec::new();
        let mut record = RecordBuf::default();
        loop {
            match reader.read_record_buf(&header, &mut record) {
                Ok(0) => break,
                Ok(_) => {
                    records.push(record.clone());
                }
                Err(e) => return Err(e.into()),
            }
        }

        Ok(records)
    }

    /// Tests that both inputs being empty produces empty output (not an error)
    #[test]
    fn test_both_inputs_empty() -> Result<()> {
        let unmapped = SamBuilder::new_unmapped();
        let mapped = SamBuilder::new_mapped();

        let records = run_zipper(&unmapped, &mapped, vec![], vec![], vec![])?;

        // Both empty should produce empty output
        assert_eq!(records.len(), 0);

        Ok(())
    }

    /// Tests the --exclude-missing-reads flag
    ///
    /// When enabled, unmapped reads without matching mapped reads should be excluded from output
    #[test]
    fn test_exclude_missing_reads() -> Result<()> {
        let mut unmapped = SamBuilder::new_unmapped();
        let mut mapped = SamBuilder::new_mapped();

        // Add unmapped pairs - q1, q2, q3
        let mut attrs = HashMap::new();
        attrs.insert("RX", BufValue::from("ACGT".to_string()));
        unmapped.add_pair_with_attrs("q1", None, None, true, true, &attrs);
        unmapped.add_pair_with_attrs("q2", None, None, true, true, &attrs);
        unmapped.add_pair_with_attrs("q3", None, None, true, true, &attrs);

        // Add mapped pairs - only q1 and q3 (q2 missing)
        let mut mapped_attrs = HashMap::new();
        mapped_attrs.insert("PG", BufValue::from(MAPPED_PG_ID.to_string()));
        mapped.add_pair_with_attrs("q1", Some(100), Some(200), true, true, &mapped_attrs);
        mapped.add_pair_with_attrs("q3", Some(300), Some(400), true, true, &mapped_attrs);

        // Without exclude_missing_reads, q2 (unmapped-only) should be included
        let records_with_missing = run_zipper_with_options(&unmapped, &mapped, false)?;
        assert_eq!(records_with_missing.len(), 6, "Should have 6 records (3 pairs)");

        // With exclude_missing_reads, q2 should be excluded
        let records_without_missing = run_zipper_with_options(&unmapped, &mapped, true)?;
        assert_eq!(records_without_missing.len(), 4, "Should have 4 records (2 pairs)");

        // Verify q2 is not in output
        for rec in &records_without_missing {
            if let Some(name) = rec.name() {
                let name_str = String::from_utf8_lossy(name.as_ref());
                assert_ne!(name_str, "q2", "q2 should be excluded");
            }
        }

        Ok(())
    }

    /// Tests that file-not-found errors are clear
    #[test]
    fn test_file_not_found_error() -> Result<()> {
        let dir = TempDir::new()?;
        let unmapped_path = dir.path().join("nonexistent_unmapped.bam");
        let mapped_path = dir.path().join("mapped.sam");
        let output_path = dir.path().join("output.bam");

        // Create mapped file but not unmapped or reference
        let mapped = SamBuilder::new_mapped();
        mapped.write_sam(&mapped_path)?;

        let zipper = Zipper {
            input: mapped_path,
            unmapped: unmapped_path.clone(),
            reference: dir.path().join("ref.fa"),
            output: output_path,
            tags_to_remove: vec![],
            tags_to_reverse: vec![],
            tags_to_revcomp: vec![],
            buffer: 5000,
            threads: 1,
            compression: CompressionOptions { compression_level: 1 },
            bwa_chunk_size: 150_000_000,
            exclude_missing_reads: false,
        };

        let result = zipper.execute("test");
        assert!(result.is_err(), "Should error when reference file doesn't exist");
        let err_msg = result.unwrap_err().to_string();
        assert!(
            err_msg.contains("does not exist") || err_msg.contains("not found"),
            "Error should mention file not found: {err_msg}"
        );

        Ok(())
    }
}

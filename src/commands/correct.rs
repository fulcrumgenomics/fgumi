//! UMI correction tool for BAM files with fixed UMI sets.
//!
//! This module provides functionality to correct Unique Molecular Identifiers (UMIs)
//! in BAM files when a fixed set of known UMI sequences is used. It can error-correct
//! observed UMIs to match the expected set based on mismatch tolerance and minimum
//! distance requirements.
//!
//! # Example
//!
//! ```no_run
//! use std::path::PathBuf;
//! use fgumi_lib::correct_umis::CorrectUmis;
//! use fgumi_lib::commands::common::{
//!     BamIoOptions, CompressionOptions, QueueMemoryOptions, RejectsOptions,
//!     SchedulerOptions, ThreadingOptions,
//! };
//!
//! let corrector = CorrectUmis {
//!     io: BamIoOptions {
//!         input: PathBuf::from("input.bam"),
//!         output: PathBuf::from("corrected.bam"),
//!     },
//!     rejects_opts: RejectsOptions { rejects: Some(PathBuf::from("rejects.bam")) },
//!     metrics: Some(PathBuf::from("metrics.txt")),
//!     max_mismatches: 2,
//!     min_distance_diff: 2,
//!     umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
//!     umi_files: vec![],
//!     umi_tag: "RX".to_string(),
//!     dont_store_original_umis: false,
//!     cache_size: 100_000,
//!     min_corrected: None,
//!     revcomp: false,
//!     threading: ThreadingOptions::new(4),
//!     compression: CompressionOptions::default(),
//!     scheduler_opts: SchedulerOptions::default(),
//!     queue_memory: QueueMemoryOptions::default(),
//! };
//!
//! corrector.execute().expect("Failed to correct UMIs");
//! ```

use ahash::AHashMap;
use anyhow::{Result, bail};
use clap::Parser;
use crossbeam_queue::SegQueue;
use fgumi_lib::bam_io::{
    create_bam_reader, create_bam_reader_for_pipeline, create_bam_writer,
    create_optional_bam_writer,
};
use fgumi_lib::bitenc::BitEnc;
use fgumi_lib::dna::reverse_complement_str;
use fgumi_lib::grouper::TemplateGrouper;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::metrics::correct::UmiCorrectionMetrics;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::template::{TemplateBatch, TemplateIterator};
use fgumi_lib::unified_pipeline::{
    BamPipelineConfig, Grouper, MemoryEstimate, run_bam_pipeline_from_reader,
    serialize_bam_records_into,
};
use fgumi_lib::validation::validate_file_exists;
use log::{info, warn};
use lru::LruCache;
use noodles::sam::alignment::io::Write as SamWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::{self, Header};
use std::io;
use std::num::NonZero;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::commands::command::Command;
use crate::commands::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, RejectsOptions, SchedulerOptions,
    ThreadingOptions,
};

/// Result of matching an observed UMI to an expected UMI.
///
/// Contains information about whether the match was acceptable and the details
/// of the best matching UMI.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UmiMatch {
    /// Whether the match is acceptable based on mismatch and distance thresholds.
    pub matched: bool,

    /// The fixed UMI sequence that was the closest match.
    pub umi: String,

    /// The number of mismatches between the observed and best matching UMI.
    pub mismatches: usize,
}

/// Corrects UMIs in BAM files to a fixed set of known UMI sequences.
///
/// This tool error-corrects Unique Molecular Identifiers (UMIs) stored in BAM files
/// when a fixed set of UMI sequences is used. It matches observed UMIs to the expected
/// set based on:
/// - Maximum allowed mismatches (`max_mismatches`)
/// - Minimum difference to the next-best match (`min_distance_diff`)
///
/// # Features
///
/// - Single and duplex (hyphenated) UMI support
/// - Optional reverse complementing of UMIs
/// - Separate output file for rejected reads
/// - Detailed metrics output in TSV format
/// - LRU caching for performance
/// - Optional storage of original UMIs before correction
/// - Parallel processing support
///
/// # Example
///
/// ```no_run
/// # use fgumi_lib::correct_umis::CorrectUmis;
/// # use fgumi_lib::commands::common::{
/// #     BamIoOptions, CompressionOptions, QueueMemoryOptions, RejectsOptions,
/// #     SchedulerOptions, ThreadingOptions,
/// # };
/// # use std::path::PathBuf;
/// let corrector = CorrectUmis {
///     io: BamIoOptions {
///         input: PathBuf::from("input.bam"),
///         output: PathBuf::from("corrected.bam"),
///     },
///     rejects_opts: RejectsOptions { rejects: Some(PathBuf::from("rejects.bam")) },
///     metrics: Some(PathBuf::from("metrics.txt")),
///     max_mismatches: 2,
///     min_distance_diff: 2,
///     umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
///     umi_files: vec![],
///     umi_tag: "RX".to_string(),
///     dont_store_original_umis: false,
///     cache_size: 100_000,
///     min_corrected: None,
///     revcomp: false,
///     threading: ThreadingOptions::new(4),
///     compression: CompressionOptions::default(),
///     scheduler_opts: SchedulerOptions::default(),
///     queue_memory: QueueMemoryOptions::default(),
/// };
///
/// corrector.execute("test")?;
/// # Ok::<(), anyhow::Error>(())
/// ```
#[derive(Parser, Debug, Clone)]
#[command(
    name = "correct",
    author,
    version,
    about = "\x1b[38;5;30m[UMI EXTRACTION]\x1b[0m \x1b[36mCorrect UMIs in a BAM file to a fixed set of UMIs\x1b[0m",
    long_about = r#"
Corrects UMIs stored in BAM files when a set of fixed UMIs is in use.

If the set of UMIs used in an experiment is known and is a _subset_ of the possible randomers
of the same length, it is possible to error-correct UMIs prior to grouping reads by UMI. This
tool takes an input BAM with UMIs in a tag (`RX` by default) and set of known UMIs (either on
the command line or in a file) and produces:

1. A new BAM with corrected UMIs in the same tag the UMIs were found in
2. Optionally a set of metrics about the representation of each UMI in the set
3. Optionally a second BAM file of reads whose UMIs could not be corrected within the specific parameters

All of the fixed UMIs must be of the same length, and all UMIs in the BAM file must also have
the same length. Multiple UMIs that are concatenated with hyphens (e.g. `AACCAGT-AGGTAGA`) are
split apart, corrected individually and then re-assembled. A read is accepted only if all the
UMIs can be corrected.

## Correction Parameters

Correction is controlled by two parameters that are applied per-UMI:

1. **--max-mismatches** controls how many mismatches (no-calls are counted as mismatches) are
   tolerated between a UMI as read and a fixed UMI
2. **--min-distance** controls how many *more* mismatches the next best hit must have

For example, with two fixed UMIs `AAAAA` and `CCCCC` and `--max-mismatches=3` and `--min-distance=2`:

- AAAAA would match to AAAAA
- AAGTG would match to AAAAA with three mismatches because CCCCC has six mismatches and 6 >= 3 + 2
- AACCA would be rejected because it is 2 mismatches to AAAAA and 3 to CCCCC and 3 <= 2 + 2

## Specifying UMIs

The set of fixed UMIs may be specified on the command line using `--umis umi1 umi2 ...` or via
one or more files of UMIs with a single sequence per line using `--umi-files umis.txt more_umis.txt`.
If there are multiple UMIs per template, leading to hyphenated UMI tags, the values for the fixed
UMIs should be single, non-hyphenated UMIs (e.g. if a record has `RX:Z:ACGT-GGCA`, you would use
`--umis ACGT GGCA`).

## Original UMI Storage

Records which have their UMIs corrected (i.e. the UMI is not identical to one of the expected
UMIs but is close enough to be corrected) will by default have their original UMI stored in the
`OX` tag. This can be disabled with the `--dont-store-original-umis` option.
"#
)]
pub struct CorrectUmis {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Options for writing rejected reads
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Optional output path for metrics TSV file.
    #[arg(short = 'm', long)]
    pub metrics: Option<PathBuf>,

    /// Maximum number of mismatches allowed.
    #[arg(long, default_value = "2")]
    pub max_mismatches: usize,

    /// Minimum difference between best and second-best match.
    #[arg(short = 'd', long = "min-distance")]
    pub min_distance_diff: usize,

    /// Fixed UMI sequences (can be specified multiple times).
    #[arg(short = 'u', long)]
    pub umis: Vec<String>,

    /// Files containing UMI sequences, one per line.
    #[arg(short = 'U', long)]
    pub umi_files: Vec<PathBuf>,

    /// SAM tag for the UMI field.
    #[arg(long, default_value = "RX")]
    pub umi_tag: String,

    /// Don't store original UMIs in a separate tag.
    #[arg(long)]
    pub dont_store_original_umis: bool,

    /// Size of the LRU cache for UMI matching.
    #[arg(long, default_value = "100000")]
    pub cache_size: usize,

    /// Minimum fraction of reads that must pass correction.
    #[arg(long)]
    pub min_corrected: Option<f64>,

    /// Reverse complement UMIs before matching.
    #[arg(long)]
    pub revcomp: bool,

    /// Threading options for parallel processing.
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline stats options
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
enum RejectionReason {
    WrongLength,
    Mismatched,
    #[default]
    None,
}

/// Result of UMI correction for a template.
///
/// This struct holds the result of correcting a UMI once for an entire template,
/// which can then be applied to all records in that template.
#[derive(Debug)]
struct TemplateCorrection {
    /// Whether the UMI matched successfully.
    matched: bool,
    /// The corrected UMI string (if matched).
    corrected_umi: Option<String>,
    /// The original UMI string.
    original_umi: String,
    /// Whether correction was needed (mismatches > 0 or revcomp).
    needs_correction: bool,
    /// Whether there were actual mismatches (not just revcomp).
    has_mismatches: bool,
    /// Match details for metrics.
    matches: Vec<UmiMatch>,
    /// Rejection reason if not matched.
    rejection_reason: RejectionReason,
}

// ============================================================================
// 7-Step Pipeline Types
// ============================================================================

/// Result from processing a batch of templates through UMI correction.
struct CorrectProcessedBatch {
    /// Records to write (kept records with corrected UMIs).
    kept_records: Vec<RecordBuf>,
    /// Rejected records (if tracking rejects).
    rejected_records: Vec<RecordBuf>,
    /// Raw-byte kept records (raw-byte mode).
    kept_raw_records: Vec<Vec<u8>>,
    /// Raw-byte rejected records (raw-byte mode).
    rejected_raw_records: Vec<Vec<u8>>,
    /// Number of templates processed.
    templates_count: u64,
    /// Number of missing UMI records.
    missing_umis: u64,
    /// Number of wrong length UMI records.
    wrong_length: u64,
    /// Number of mismatched UMI records.
    mismatched: u64,
    /// Per-UMI match counts for metrics.
    umi_matches: AHashMap<String, UmiCorrectionMetrics>,
}

impl MemoryEstimate for CorrectProcessedBatch {
    fn estimate_heap_size(&self) -> usize {
        let recordbuf_size: usize = self
            .kept_records
            .iter()
            .chain(self.rejected_records.iter())
            .map(MemoryEstimate::estimate_heap_size)
            .sum();
        let recordbuf_vec_overhead = (self.kept_records.capacity()
            + self.rejected_records.capacity())
            * std::mem::size_of::<RecordBuf>();
        let raw_size: usize = self
            .kept_raw_records
            .iter()
            .chain(self.rejected_raw_records.iter())
            .map(Vec::capacity)
            .sum();
        let raw_vec_overhead = (self.kept_raw_records.capacity()
            + self.rejected_raw_records.capacity())
            * std::mem::size_of::<Vec<u8>>();
        recordbuf_size + recordbuf_vec_overhead + raw_size + raw_vec_overhead
    }
}

/// Metrics collected from UMI correction processing, aggregated post-pipeline.
#[derive(Default)]
struct CollectedCorrectMetrics {
    /// Total templates processed.
    templates_processed: u64,
    /// Records with missing UMI tag.
    missing_umis: u64,
    /// Records with wrong UMI length.
    wrong_length: u64,
    /// Records that didn't match any fixed UMI.
    mismatched: u64,
    /// Per-UMI match counts (for metrics file).
    umi_matches: AHashMap<String, UmiCorrectionMetrics>,
    /// Rejected records (if tracking).
    rejects: Vec<RecordBuf>,
    /// Rejected raw records (if tracking, for raw-byte mode).
    raw_rejects: Vec<Vec<u8>>,
}

impl Command for CorrectUmis {
    /// Executes the UMI correction pipeline.
    ///
    /// This method:
    /// 1. Validates input parameters
    /// 2. Loads UMI sequences
    /// 3. Checks UMI distances and warns about ambiguities
    /// 4. Processes the BAM file in parallel
    /// 5. Writes corrected/rejected records
    /// 6. Computes and outputs metrics
    ///
    /// # Returns
    ///
    /// `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Validation fails
    /// - Files cannot be read/written
    /// - BAM processing encounters errors
    /// - Metrics cannot be written
    /// - Minimum correction ratio is not met
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use fgumi_lib::correct_umis::CorrectUmis;
    /// # use fgumi_lib::commands::common::{BamIoOptions, RejectsOptions, ThreadingOptions};
    /// # use std::path::PathBuf;
    /// # use fgumi_lib::commands::command::Command;
    /// let corrector = CorrectUmis {
    ///     /* ... field initialization ... */
    /// #   io: BamIoOptions { input: PathBuf::new(), output: PathBuf::new() },
    /// #   rejects_opts: RejectsOptions::default(),
    /// #   metrics: None,
    /// #   max_mismatches: 2,
    /// #   min_distance_diff: 2,
    /// #   umis: vec![],
    /// #   umi_files: vec![],
    /// #   umi_tag: "RX".to_string(),
    /// #   dont_store_original_umis: false,
    /// #   cache_size: 100_000,
    /// #   min_corrected: None,
    /// #   revcomp: false,
    /// #   threading: ThreadingOptions::new(4),
    /// };
    ///
    /// corrector.execute("test")?;
    /// # Ok::<(), anyhow::Error>(())
    /// ```
    fn execute(&self, command_line: &str) -> Result<()> {
        self.validate()?;

        let timer = OperationTimer::new("Correcting UMIs");

        // Load UMI sequences
        let (umi_sequences, umi_length) = self.load_umi_sequences()?;
        // Create encoded UMI set for fast comparison
        let encoded_umi_set = EncodedUmiSet::new(&umi_sequences);

        // Warn about UMIs that are too close together
        self.check_umi_distances(&umi_sequences);

        // Open input using streaming-capable reader for pipeline use
        let (reader, header) = create_bam_reader_for_pipeline(&self.io.input)?;

        // Add @PG record with PP chaining to input's last program
        let header = fgumi_lib::header::add_pg_record(
            header,
            crate::version::VERSION.as_str(),
            command_line,
        )?;

        // Dispatch based on thread count:
        // - Some(threads) -> 7-step pipeline
        // - None -> single-threaded fast path
        let total_records = if let Some(threads) = self.threading.threads {
            // Multi-threaded: 7-step pipeline
            self.execute_threads_mode(
                threads,
                reader,
                header,
                Arc::new(encoded_umi_set),
                umi_length,
                self.rejects_opts.rejects.is_some(),
            )?
        } else {
            // Drop reader for single-threaded mode - uses its own iterator
            drop(reader);
            // Single-threaded: main thread only, template-level optimization
            self.execute_single_thread_mode(
                header,
                encoded_umi_set,
                umi_length,
                self.rejects_opts.rejects.is_some(),
            )?
        };

        timer.log_completion(total_records);
        Ok(())
    }
}

impl CorrectUmis {
    /// Validates input parameters.
    ///
    /// Checks that:
    /// - At least one UMI is provided
    /// - Input file exists
    /// - `min_corrected` is between 0 and 1 if specified
    ///
    /// # Errors
    ///
    /// Returns an error if any validation check fails.
    fn validate(&self) -> Result<()> {
        if self.umis.is_empty() && self.umi_files.is_empty() {
            bail!("At least one UMI or UMI file must be provided.");
        }
        validate_file_exists(&self.io.input, "input BAM file")?;
        if let Some(min) = self.min_corrected {
            if !(0.0..=1.0).contains(&min) {
                bail!("--min-corrected must be between 0 and 1.");
            }
        }
        Ok(())
    }

    /// Loads UMI sequences from command line and files.
    ///
    /// Combines UMIs from both `umis` and `umi_files`, converts them to uppercase,
    /// and validates that all UMIs are the same length.
    ///
    /// # Returns
    ///
    /// A tuple of `(umi_sequences, umi_length)`.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - No UMIs are provided
    /// - UMI files cannot be read
    /// - UMIs have different lengths
    fn load_umi_sequences(&self) -> Result<(Vec<String>, usize)> {
        let mut umi_set: std::collections::HashSet<String> =
            self.umis.iter().map(|s| s.to_uppercase()).collect();

        for file in &self.umi_files {
            let content = std::fs::read_to_string(file)?;
            for line in content.lines() {
                let umi = line.trim().to_uppercase();
                if !umi.is_empty() {
                    umi_set.insert(umi);
                }
            }
        }

        if umi_set.is_empty() {
            bail!("No UMIs provided.");
        }

        let mut umi_sequences: Vec<String> = umi_set.into_iter().collect();
        umi_sequences.sort_unstable();

        // Check all UMIs have the same length
        let first_len = umi_sequences[0].len();
        if !umi_sequences.iter().all(|u| u.len() == first_len) {
            bail!("All UMIs must have the same length.");
        }

        info!("Loaded {} UMI sequences of length {}", umi_sequences.len(), first_len);
        Ok((umi_sequences, first_len))
    }

    /// Checks distances between UMI pairs and warns about ambiguities.
    ///
    /// Identifies pairs of UMIs that are within `min_distance_diff` of each other,
    /// which could lead to ambiguous matching situations.
    ///
    /// # Arguments
    ///
    /// * `umi_sequences` - Slice of UMI sequences to check
    fn check_umi_distances(&self, umi_sequences: &[String]) {
        let pairs = find_umi_pairs_within_distance(umi_sequences, self.min_distance_diff - 1);

        if !pairs.is_empty() {
            warn!("###################################################################");
            warn!("# WARNING: Found pairs of UMIs within min-distance-diff threshold!");
            warn!("# These pairs may be ambiguous and fail to match:");
            for (umi1, umi2, dist) in &pairs {
                warn!("#   {umi1} <-> {umi2} (distance {dist})");
            }
            warn!("###################################################################");
        }
    }

    /// Extracts and validates that all records in a template have the same UMI.
    ///
    /// Returns:
    /// - `Some(umi)` if all records have the same UMI
    /// - `None` if no records have a UMI (all missing)
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - Records have different UMIs
    /// - Some records have UMIs and others don't
    /// - UMI tag has non-string value
    fn extract_and_validate_template_umi(records: &[RecordBuf], umi_tag: Tag) -> Option<String> {
        if records.is_empty() {
            return None;
        }

        // Extract UMI from first record
        let first_umi = records[0].data().get(&umi_tag).map(|v| {
            if let sam::alignment::record_buf::data::field::Value::String(s) = v {
                String::from_utf8_lossy(s).to_string()
            } else {
                panic!(
                    "UMI tag {:?} has non-string value in record {:?}",
                    umi_tag,
                    records[0].name()
                );
            }
        });

        // Validate all other records have the same UMI
        for record in &records[1..] {
            let current_umi = record.data().get(&umi_tag).map(|v| {
                if let sam::alignment::record_buf::data::field::Value::String(s) = v {
                    String::from_utf8_lossy(s).to_string()
                } else {
                    panic!(
                        "UMI tag {:?} has non-string value in record {:?}",
                        umi_tag,
                        record.name()
                    );
                }
            });

            match (&first_umi, &current_umi) {
                (Some(first), Some(current)) if first != current => {
                    panic!(
                        "Template {:?} has mismatched UMIs: first={:?}, current={:?}",
                        records[0].name(),
                        first,
                        current
                    );
                }
                (Some(_), None) | (None, Some(_)) => {
                    panic!(
                        "Template {:?} has inconsistent UMI presence across records",
                        records[0].name()
                    );
                }
                _ => {} // Both None or both Some with matching values
            }
        }

        first_umi
    }

    /// Compute UMI correction for a template (called once per template).
    #[allow(clippy::too_many_arguments)]
    fn compute_template_correction(
        umi: &str,
        umi_length: usize,
        revcomp: bool,
        max_mismatches: usize,
        min_distance_diff: usize,
        encoded_umi_set: &EncodedUmiSet,
        cache: &mut Option<LruCache<Vec<u8>, UmiMatch>>,
    ) -> TemplateCorrection {
        let original_umi = umi.to_string();

        // Split and optionally reverse complement
        let sequences: Vec<String> = if revcomp {
            umi.split('-').map(reverse_complement_str).rev().collect()
        } else {
            umi.split('-').map(std::string::ToString::to_string).collect()
        };

        // Check length
        if sequences.iter().any(|s| s.len() != umi_length) {
            return TemplateCorrection {
                matched: false,
                corrected_umi: None,
                original_umi,
                needs_correction: false,
                has_mismatches: false,
                matches: Vec::new(),
                rejection_reason: RejectionReason::WrongLength,
            };
        }

        // Find matches for each UMI segment
        let mut matches = Vec::with_capacity(sequences.len());
        for seq in &sequences {
            // Uppercase using bytes directly - avoids String allocation
            let seq_bytes: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();

            let umi_match = if let Some(c) = cache {
                if let Some(cached) = c.get(&seq_bytes[..]) {
                    cached.clone()
                } else {
                    let result = find_best_match_encoded(
                        &seq_bytes,
                        encoded_umi_set,
                        max_mismatches,
                        min_distance_diff,
                    );
                    c.put(seq_bytes, result.clone());
                    result
                }
            } else {
                find_best_match_encoded(
                    &seq_bytes,
                    encoded_umi_set,
                    max_mismatches,
                    min_distance_diff,
                )
            };

            matches.push(umi_match);
        }

        // Determine if all segments matched
        let all_matched = matches.iter().all(|m| m.matched);
        let has_mismatches = matches.iter().any(|m| m.mismatches > 0);
        let needs_correction = has_mismatches || revcomp;

        if all_matched {
            let corrected_umi: String =
                matches.iter().map(|m| m.umi.clone()).collect::<Vec<_>>().join("-");
            TemplateCorrection {
                matched: true,
                corrected_umi: Some(corrected_umi),
                original_umi,
                needs_correction,
                has_mismatches,
                matches,
                rejection_reason: RejectionReason::None,
            }
        } else {
            TemplateCorrection {
                matched: false,
                corrected_umi: None,
                original_umi,
                needs_correction: false,
                has_mismatches: false,
                matches,
                rejection_reason: RejectionReason::Mismatched,
            }
        }
    }

    /// Apply a computed correction to a record.
    fn apply_correction_to_record(
        mut record: RecordBuf,
        correction: &TemplateCorrection,
        umi_tag: Tag,
        dont_store_original_umis: bool,
    ) -> RecordBuf {
        if correction.needs_correction {
            // Store original UMI if there were actual mismatches
            if !dont_store_original_umis && correction.has_mismatches {
                record.data_mut().insert(
                    Tag::ORIGINAL_UMI_BARCODE_SEQUENCE,
                    sam::alignment::record_buf::data::field::Value::String(
                        correction.original_umi.clone().into(),
                    ),
                );
            }

            // Write corrected UMI
            if let Some(ref corrected) = correction.corrected_umi {
                record.data_mut().insert(
                    umi_tag,
                    sam::alignment::record_buf::data::field::Value::String(
                        corrected.clone().into(),
                    ),
                );
            }
        }
        record
    }

    /// Extract and validate UMI from raw-byte records in a template.
    fn extract_and_validate_template_umi_raw(
        raw_records: &[Vec<u8>],
        umi_tag: [u8; 2],
    ) -> Option<String> {
        use fgumi_lib::sort::bam_fields;

        if raw_records.is_empty() {
            return None;
        }

        // Guard against truncated records
        if raw_records.iter().any(|r| r.len() < 32) {
            return None;
        }

        let first_aux = bam_fields::aux_data_slice(&raw_records[0]);
        // Work with &[u8] slices to avoid per-record String allocation
        let first_umi_bytes = bam_fields::find_string_tag(first_aux, &umi_tag);

        // If tag exists but is not a Z-type string, panic (matching RecordBuf behavior)
        if first_umi_bytes.is_none() {
            if let Some(tag_type) = bam_fields::find_tag_type(first_aux, &umi_tag) {
                panic!(
                    "UMI tag {:?} exists but has non-string type '{}', expected 'Z'",
                    std::str::from_utf8(&umi_tag).unwrap_or("??"),
                    tag_type as char,
                );
            }
        }

        for raw in &raw_records[1..] {
            let aux = bam_fields::aux_data_slice(raw);
            let current_umi_bytes = bam_fields::find_string_tag(aux, &umi_tag);

            match (first_umi_bytes, current_umi_bytes) {
                (Some(first), Some(current)) if first != current => {
                    panic!(
                        "Template has mismatched UMIs: first={:?}, current={:?}",
                        String::from_utf8_lossy(first),
                        String::from_utf8_lossy(current)
                    );
                }
                (Some(_), None) | (None, Some(_)) => {
                    panic!("Template has inconsistent UMI presence across records");
                }
                _ => {}
            }
        }

        // Only allocate String once at the end
        first_umi_bytes.map(|b| String::from_utf8_lossy(b).into_owned())
    }

    /// Apply UMI correction to a raw BAM record.
    fn apply_correction_to_raw(
        record: &mut Vec<u8>,
        correction: &TemplateCorrection,
        umi_tag: [u8; 2],
        dont_store_original_umis: bool,
    ) {
        use fgumi_lib::sort::bam_fields;

        if correction.needs_correction {
            // Write corrected UMI first (in-place update avoids scanning past OX)
            if let Some(ref corrected) = correction.corrected_umi {
                bam_fields::update_string_tag(record, &umi_tag, corrected.as_bytes());
            }

            // Store original UMI if there were actual mismatches
            // Use update_string_tag to avoid duplicate OX tags
            if !dont_store_original_umis && correction.has_mismatches {
                bam_fields::update_string_tag(record, b"OX", correction.original_umi.as_bytes());
            }
        }
    }

    fn finalize_metrics(
        &self,
        umi_metrics: &mut AHashMap<String, UmiCorrectionMetrics>,
        unmatched_umi: &str,
    ) -> Result<()> {
        // Calculate totals
        let total: u64 = umi_metrics.values().map(|m| m.total_matches).sum();
        let matched_total: u64 = umi_metrics
            .iter()
            .filter(|(umi, _)| *umi != unmatched_umi)
            .map(|(_, m)| m.total_matches)
            .sum();

        // Calculate fractions (allow NaN when total is 0, matching fgbio behavior)
        #[allow(clippy::cast_precision_loss)]
        for metric in umi_metrics.values_mut() {
            metric.fraction_of_matches = metric.total_matches as f64 / total as f64;
        }

        // Calculate representations (allow NaN/Infinity, matching fgbio behavior)
        let umi_count = umi_metrics.keys().filter(|umi| *umi != unmatched_umi).count();

        #[allow(clippy::cast_precision_loss)]
        let mean = matched_total as f64 / umi_count as f64;
        for metric in umi_metrics.values_mut() {
            metric.representation = metric.total_matches as f64 / mean;
        }

        // Write metrics if requested
        if let Some(path) = &self.metrics {
            let mut metrics: Vec<UmiCorrectionMetrics> = umi_metrics.values().cloned().collect();
            metrics.sort_by(|a, b| a.umi.cmp(&b.umi));
            UmiCorrectionMetrics::write_metrics(&metrics, path)?;
        }

        Ok(())
    }

    /// Execute using the 7-step unified pipeline.
    ///
    /// This method uses the template-based grouper to process records in batches,
    /// with parallel decompression, processing, and compression.
    #[allow(clippy::too_many_lines)]
    fn execute_threads_mode(
        &self,
        num_threads: usize,
        reader: Box<dyn std::io::Read + Send>,
        header: Header,
        encoded_umi_set: Arc<EncodedUmiSet>,
        umi_length: usize,
        track_rejects: bool,
    ) -> Result<u64> {
        // Configure pipeline - correct is ReaderHeavy (70% in decompression)
        let mut pipeline_config =
            BamPipelineConfig::auto_tuned(num_threads, self.compression.compression_level);
        pipeline_config.pipeline.scheduler_strategy = self.scheduler_opts.strategy();
        if self.scheduler_opts.collect_stats() {
            pipeline_config.pipeline = pipeline_config.pipeline.with_stats(true);
        }
        pipeline_config.pipeline.deadlock_timeout_secs =
            self.scheduler_opts.deadlock_timeout_secs();
        pipeline_config.pipeline.deadlock_recover_enabled =
            self.scheduler_opts.deadlock_recover_enabled();

        // Calculate and apply queue memory limit
        let queue_memory_limit_bytes = self.queue_memory.calculate_memory_limit(num_threads)?;
        pipeline_config.pipeline.queue_memory_limit = queue_memory_limit_bytes;
        self.queue_memory.log_memory_config(num_threads, queue_memory_limit_bytes);

        // Enable raw-byte mode (correct uses TemplateGrouper, no cell tag needed)
        {
            use fgumi_lib::read_info::LibraryIndex;
            use fgumi_lib::unified_pipeline::GroupKeyConfig;

            let library_index = LibraryIndex::from_header(&header);
            pipeline_config.group_key_config = Some(GroupKeyConfig::new_raw_no_cell(library_index));
        }

        // Lock-free metrics collection
        let collected_metrics: Arc<SegQueue<CollectedCorrectMetrics>> = Arc::new(SegQueue::new());
        let collected_for_serialize = Arc::clone(&collected_metrics);

        // Configuration for closures
        const BATCH_SIZE: usize = 1000; // Templates per batch
        let max_mismatches = self.max_mismatches;
        let min_distance_diff = self.min_distance_diff;
        let umi_tag = Tag::new(self.umi_tag.as_bytes()[0], self.umi_tag.as_bytes()[1]);
        let revcomp = self.revcomp;
        let cache_size = self.cache_size;
        let dont_store_original_umis = self.dont_store_original_umis;

        // Progress tracking
        let progress_counter = Arc::new(AtomicU64::new(0));
        let progress_for_process = Arc::clone(&progress_counter);

        // Unmatched UMI string for rejected reads (matches fgbio behavior)
        let unmatched_umi = "N".repeat(umi_length);
        let unmatched_umi_for_process = unmatched_umi.clone();

        // Clone the UMI set for post-aggregation use (before closure moves original)
        let encoded_umi_set_for_metrics = Arc::clone(&encoded_umi_set);

        // Grouper: batch templates by QNAME
        let grouper_fn = move |_header: &Header| {
            Box::new(TemplateGrouper::new(BATCH_SIZE))
                as Box<dyn Grouper<Group = TemplateBatch> + Send>
        };

        // Process function: correct UMIs in each template batch
        let process_fn = move |batch: TemplateBatch| -> io::Result<CorrectProcessedBatch> {
            // Per-thread LRU cache for UMI matching
            thread_local! {
                static CACHE: std::cell::RefCell<Option<LruCache<Vec<u8>, UmiMatch>>> = const { std::cell::RefCell::new(None) };
            }

            CACHE.with(|cache_cell| {
                let mut cache_ref = cache_cell.borrow_mut();
                if cache_ref.is_none() && cache_size > 0 {
                    *cache_ref = Some(LruCache::new(NonZero::new(cache_size).unwrap()));
                }

                let mut kept_records = Vec::new();
                let mut rejected_records = Vec::new();
                let mut kept_raw_records = Vec::new();
                let mut rejected_raw_records: Vec<Vec<u8>> = Vec::new();
                let mut missing_umis = 0u64;
                let mut wrong_length = 0u64;
                let mut mismatched = 0u64;
                let mut umi_matches_map: AHashMap<String, UmiCorrectionMetrics> = AHashMap::new();
                let templates_count = batch.len() as u64;
                // Count ALL input records for progress tracking (not just kept/rejected)
                let mut total_input_records = 0u64;
                let umi_tag_bytes: [u8; 2] = [umi_tag.as_ref()[0], umi_tag.as_ref()[1]];

                for template in batch {
                    // Count input records BEFORE processing
                    total_input_records += template.read_count() as u64;

                    if template.is_raw_byte_mode() {
                        // Raw-byte mode - take ownership to avoid deep clone
                        let raw_records = template.into_raw_records().unwrap_or_default();
                        let umi_opt = Self::extract_and_validate_template_umi_raw(
                            &raw_records,
                            umi_tag_bytes,
                        );

                        match umi_opt {
                            None => {
                                let num_records = raw_records.len() as u64;
                                missing_umis += num_records;
                                let entry = umi_matches_map
                                    .entry(unmatched_umi_for_process.clone())
                                    .or_insert_with(|| {
                                        UmiCorrectionMetrics::new(unmatched_umi_for_process.clone())
                                    });
                                entry.total_matches += num_records;
                                if track_rejects {
                                    rejected_raw_records.extend(raw_records);
                                }
                            }
                            Some(umi) => {
                                let correction = Self::compute_template_correction(
                                    &umi,
                                    umi_length,
                                    revcomp,
                                    max_mismatches,
                                    min_distance_diff,
                                    &encoded_umi_set,
                                    &mut cache_ref,
                                );
                                let num_records = raw_records.len() as u64;

                                if correction.matched {
                                    for m in &correction.matches {
                                        if m.matched {
                                            let entry = umi_matches_map
                                                .entry(m.umi.clone())
                                                .or_insert_with(|| {
                                                    UmiCorrectionMetrics::new(m.umi.clone())
                                                });
                                            entry.total_matches += num_records;
                                            match m.mismatches {
                                                0 => entry.perfect_matches += num_records,
                                                1 => entry.one_mismatch_matches += num_records,
                                                2 => entry.two_mismatch_matches += num_records,
                                                _ => entry.other_matches += num_records,
                                            }
                                        }
                                    }

                                    for mut raw in raw_records {
                                        Self::apply_correction_to_raw(
                                            &mut raw,
                                            &correction,
                                            umi_tag_bytes,
                                            dont_store_original_umis,
                                        );
                                        kept_raw_records.push(raw);
                                    }
                                } else {
                                    match correction.rejection_reason {
                                        RejectionReason::WrongLength => {
                                            wrong_length += num_records;
                                        }
                                        RejectionReason::Mismatched => {
                                            mismatched += num_records;
                                        }
                                        RejectionReason::None => {}
                                    }
                                    let entry = umi_matches_map
                                        .entry(unmatched_umi_for_process.clone())
                                        .or_insert_with(|| {
                                            UmiCorrectionMetrics::new(
                                                unmatched_umi_for_process.clone(),
                                            )
                                        });
                                    entry.total_matches += num_records;
                                    if track_rejects {
                                        rejected_raw_records.extend(raw_records);
                                    }
                                }
                            }
                        }
                    } else {
                        // RecordBuf mode
                        let records: Vec<RecordBuf> = template.into_records();
                        let umi_opt = Self::extract_and_validate_template_umi(&records, umi_tag);

                        match umi_opt {
                            None => {
                                let num_records = records.len() as u64;
                                missing_umis += num_records;
                                let entry = umi_matches_map
                                    .entry(unmatched_umi_for_process.clone())
                                    .or_insert_with(|| {
                                        UmiCorrectionMetrics::new(unmatched_umi_for_process.clone())
                                    });
                                entry.total_matches += num_records;
                                if track_rejects {
                                    rejected_records.extend(records);
                                }
                            }
                            Some(umi) => {
                                let correction = Self::compute_template_correction(
                                    &umi,
                                    umi_length,
                                    revcomp,
                                    max_mismatches,
                                    min_distance_diff,
                                    &encoded_umi_set,
                                    &mut cache_ref,
                                );
                                let num_records = records.len() as u64;

                                if correction.matched {
                                    for m in &correction.matches {
                                        if m.matched {
                                            let entry = umi_matches_map
                                                .entry(m.umi.clone())
                                                .or_insert_with(|| {
                                                    UmiCorrectionMetrics::new(m.umi.clone())
                                                });
                                            entry.total_matches += num_records;
                                            match m.mismatches {
                                                0 => entry.perfect_matches += num_records,
                                                1 => entry.one_mismatch_matches += num_records,
                                                2 => entry.two_mismatch_matches += num_records,
                                                _ => entry.other_matches += num_records,
                                            }
                                        }
                                    }

                                    for record in records {
                                        let corrected = Self::apply_correction_to_record(
                                            record,
                                            &correction,
                                            umi_tag,
                                            dont_store_original_umis,
                                        );
                                        kept_records.push(corrected);
                                    }
                                } else {
                                    match correction.rejection_reason {
                                        RejectionReason::WrongLength => {
                                            wrong_length += num_records;
                                        }
                                        RejectionReason::Mismatched => {
                                            mismatched += num_records;
                                        }
                                        RejectionReason::None => {}
                                    }
                                    let entry = umi_matches_map
                                        .entry(unmatched_umi_for_process.clone())
                                        .or_insert_with(|| {
                                            UmiCorrectionMetrics::new(
                                                unmatched_umi_for_process.clone(),
                                            )
                                        });
                                    entry.total_matches += num_records;
                                    if track_rejects {
                                        rejected_records.extend(records);
                                    }
                                }
                            }
                        }
                    }
                }

                // Progress logging (count ALL input records, not just output)
                let count = progress_for_process.fetch_add(total_input_records, Ordering::Relaxed);
                if (count + total_input_records) / 1_000_000 > count / 1_000_000 {
                    info!("Processed {} records", count + total_input_records);
                }

                Ok(CorrectProcessedBatch {
                    kept_records,
                    rejected_records,
                    kept_raw_records,
                    rejected_raw_records,
                    templates_count,
                    missing_umis,
                    wrong_length,
                    mismatched,
                    umi_matches: umi_matches_map,
                })
            })
        };

        // Serialize function: convert records to bytes and collect metrics
        let serialize_fn = move |processed: CorrectProcessedBatch,
                                 header: &Header,
                                 output: &mut Vec<u8>|
              -> io::Result<u64> {
            // Push metrics to lock-free queue
            collected_for_serialize.push(CollectedCorrectMetrics {
                templates_processed: processed.templates_count,
                missing_umis: processed.missing_umis,
                wrong_length: processed.wrong_length,
                mismatched: processed.mismatched,
                umi_matches: processed.umi_matches,
                rejects: processed.rejected_records,
                raw_rejects: processed.rejected_raw_records,
            });

            // Serialize kept records: raw-byte path or RecordBuf path
            let mut kept_count = 0u64;
            if !processed.kept_raw_records.is_empty() {
                for raw in &processed.kept_raw_records {
                    let block_size = raw.len() as u32;
                    output.extend_from_slice(&block_size.to_le_bytes());
                    output.extend_from_slice(raw);
                    kept_count += 1;
                }
            }
            if !processed.kept_records.is_empty() {
                kept_count += processed.kept_records.len() as u64;
                serialize_bam_records_into(&processed.kept_records, header, output)?;
            }
            // Return KEPT record count (not input count, to avoid double-counting)
            Ok(kept_count)
        };

        // Clone header before pipeline (needed for post-pipeline rejects writing)
        let header_for_rejects = header.clone();

        // Run the 7-step pipeline with the already-opened reader (supports streaming)
        let records_written = run_bam_pipeline_from_reader(
            pipeline_config,
            reader,
            header,
            &self.io.output,
            None, // Use input header for output
            grouper_fn,
            process_fn,
            serialize_fn,
        )?;

        // ========== Post-pipeline: Aggregate metrics ==========
        let mut total_templates = 0u64;
        let mut total_missing = 0u64;
        let mut total_wrong_length = 0u64;
        let mut total_mismatched = 0u64;
        let mut merged_umi_matches: AHashMap<String, UmiCorrectionMetrics> = AHashMap::new();
        let mut all_rejects: Vec<RecordBuf> = Vec::new();
        let mut all_raw_rejects: Vec<Vec<u8>> = Vec::new();

        while let Some(metrics) = collected_metrics.pop() {
            total_templates += metrics.templates_processed;
            total_missing += metrics.missing_umis;
            total_wrong_length += metrics.wrong_length;
            total_mismatched += metrics.mismatched;

            for (umi, counts) in metrics.umi_matches {
                let entry = merged_umi_matches
                    .entry(umi.clone())
                    .or_insert_with(|| UmiCorrectionMetrics::new(umi));
                entry.total_matches += counts.total_matches;
                entry.perfect_matches += counts.perfect_matches;
                entry.one_mismatch_matches += counts.one_mismatch_matches;
                entry.two_mismatch_matches += counts.two_mismatch_matches;
                entry.other_matches += counts.other_matches;
            }

            if track_rejects {
                all_rejects.extend(metrics.rejects);
                all_raw_rejects.extend(metrics.raw_rejects);
            }
        }

        // Write rejects file if requested
        let has_rejects = !all_rejects.is_empty() || !all_raw_rejects.is_empty();
        if track_rejects && has_rejects {
            if let Some(rejects_path) = &self.rejects_opts.rejects {
                let writer_threads = self.threading.num_threads();
                let mut rejects_writer = create_optional_bam_writer(
                    Some(rejects_path),
                    &header_for_rejects,
                    writer_threads,
                    self.compression.compression_level,
                )?
                .expect("rejects path is Some");

                // Ensure we never mix RecordBuf and raw rejects in the same run
                debug_assert!(
                    all_rejects.is_empty() || all_raw_rejects.is_empty(),
                    "mixed RecordBuf and raw-byte rejects in same run"
                );

                // Write RecordBuf rejects
                for record in &all_rejects {
                    rejects_writer.write_alignment_record(&header_for_rejects, record)?;
                }
                // Write raw-byte rejects
                for raw in &all_raw_rejects {
                    use std::io::Write;
                    let block_size = raw.len() as u32;
                    rejects_writer.get_mut().write_all(&block_size.to_le_bytes())?;
                    rejects_writer.get_mut().write_all(raw)?;
                }

                rejects_writer.finish(&header_for_rejects)?;
                let total_reject_records = all_rejects.len() + all_raw_rejects.len();
                info!("Wrote {total_reject_records} rejected records to rejects file");
            }
        }

        // Ensure ALL UMI rows are present (even with zero counts) - matches fgbio behavior
        for umi in encoded_umi_set_for_metrics.strings.iter().chain(std::iter::once(&unmatched_umi))
        {
            merged_umi_matches
                .entry(umi.clone())
                .or_insert_with(|| UmiCorrectionMetrics::new(umi.clone()));
        }

        // Finalize and write metrics if requested
        self.finalize_metrics(&mut merged_umi_matches, &unmatched_umi)?;

        // Log summary (records_written = kept records only)
        let rejected = total_missing + total_wrong_length + total_mismatched;
        let total_records = records_written + rejected;
        info!("Read {total_records}; kept {records_written} and rejected {rejected}");
        info!("Total templates processed: {total_templates}");

        if total_missing > 0 || total_wrong_length > 0 {
            warn!("###################################################################");
            if total_missing > 0 {
                warn!("# {total_missing} were missing UMI attributes in the BAM file!");
            }
            if total_wrong_length > 0 {
                warn!(
                    "# {total_wrong_length} had unexpected UMIs of differing lengths in the BAM file!"
                );
            }
            warn!("###################################################################");
        }

        // Check minimum correction ratio
        if let Some(min) = self.min_corrected {
            #[allow(clippy::cast_precision_loss)]
            let ratio_kept = records_written as f64 / total_records as f64;
            if ratio_kept < min {
                bail!(
                    "Final ratio of reads kept / total was {ratio_kept:.2} (user specified minimum was {min:.2}). \
                    This could indicate a mismatch between library preparation and the provided UMI file."
                );
            }
        }

        Ok(total_records)
    }

    /// Execute using single-threaded mode with template-level UMI correction.
    ///
    /// This method runs entirely on the main thread with no pipeline overhead.
    /// It uses `TemplateBatchReader` to read templates and applies UMI correction
    /// once per template, then applies the correction to all records in the template.
    #[allow(clippy::too_many_lines)]
    fn execute_single_thread_mode(
        &self,
        header: Header,
        encoded_umi_set: EncodedUmiSet,
        umi_length: usize,
        track_rejects: bool,
    ) -> Result<u64> {
        info!("Using single-threaded mode with template-level UMI correction");

        // Open input BAM reader and create template iterator
        let (mut bam_reader, _) = create_bam_reader(&self.io.input, 1)?;
        let record_iter = bam_reader.record_bufs(&header).map(|r| r.map_err(Into::into));
        let template_iter = TemplateIterator::new(record_iter);

        // Open output writer (single-threaded)
        let mut writer =
            create_bam_writer(&self.io.output, &header, 1, self.compression.compression_level)?;
        let mut reject_writer = create_optional_bam_writer(
            self.rejects_opts.rejects.as_ref(),
            &header,
            1,
            self.compression.compression_level,
        )?;

        // LRU cache (single thread, no thread_local needed)
        let mut cache: Option<LruCache<Vec<u8>, UmiMatch>> = if self.cache_size > 0 {
            Some(LruCache::new(NonZero::new(self.cache_size).unwrap()))
        } else {
            None
        };

        // Initialize metrics
        let unmatched_umi = "N".repeat(umi_length);
        let umi_sequences = encoded_umi_set.strings.clone();
        let mut umi_metrics: AHashMap<String, UmiCorrectionMetrics> = umi_sequences
            .iter()
            .chain(std::iter::once(&unmatched_umi))
            .map(|umi| (umi.clone(), UmiCorrectionMetrics::new(umi.clone())))
            .collect();

        // Counters
        let mut total_records = 0u64;
        let mut total_templates = 0u64;
        let mut missing_umis = 0u64;
        let mut wrong_length = 0u64;
        let mut mismatched = 0u64;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        let umi_tag = Tag::new(self.umi_tag.as_bytes()[0], self.umi_tag.as_bytes()[1]);
        let max_mismatches = self.max_mismatches;
        let min_distance_diff = self.min_distance_diff;
        let revcomp = self.revcomp;
        let dont_store_original_umis = self.dont_store_original_umis;

        // Process templates
        for template_result in template_iter {
            let template = template_result?;
            let records: Vec<RecordBuf> = template.into_records();
            let num_records = records.len() as u64;
            total_records += num_records;
            total_templates += 1;

            // Template-level UMI correction
            let umi_opt = Self::extract_and_validate_template_umi(&records, umi_tag);

            match umi_opt {
                None => {
                    // No UMI in any record - reject all records
                    missing_umis += num_records;
                    umi_metrics.get_mut(&unmatched_umi).unwrap().total_matches += num_records;
                    if track_rejects {
                        if let Some(ref mut rw) = reject_writer {
                            for record in records {
                                rw.write_alignment_record(&header, &record)?;
                            }
                        }
                    }
                }
                Some(umi) => {
                    // Correct UMI once for the template
                    let correction = Self::compute_template_correction(
                        &umi,
                        umi_length,
                        revcomp,
                        max_mismatches,
                        min_distance_diff,
                        &encoded_umi_set,
                        &mut cache,
                    );

                    if correction.matched {
                        // Update metrics for matched UMIs (count per-read, not per-template)
                        for m in &correction.matches {
                            if m.matched {
                                let entry = umi_metrics.get_mut(&m.umi).unwrap();
                                entry.total_matches += num_records;
                                match m.mismatches {
                                    0 => entry.perfect_matches += num_records,
                                    1 => entry.one_mismatch_matches += num_records,
                                    2 => entry.two_mismatch_matches += num_records,
                                    _ => entry.other_matches += num_records,
                                }
                            }
                        }

                        // Apply correction to all records in template and write
                        for record in records {
                            let corrected = Self::apply_correction_to_record(
                                record,
                                &correction,
                                umi_tag,
                                dont_store_original_umis,
                            );
                            writer.write_alignment_record(&header, &corrected)?;
                        }
                    } else {
                        // Rejection - update counters based on reason
                        match correction.rejection_reason {
                            RejectionReason::WrongLength => wrong_length += num_records,
                            RejectionReason::Mismatched => mismatched += num_records,
                            RejectionReason::None => {}
                        }
                        umi_metrics.get_mut(&unmatched_umi).unwrap().total_matches += num_records;

                        if track_rejects {
                            if let Some(ref mut rw) = reject_writer {
                                for record in records {
                                    rw.write_alignment_record(&header, &record)?;
                                }
                            }
                        }
                    }
                }
            }

            progress.log_if_needed(num_records);
        }

        progress.log_final();

        // Finish writers
        writer.finish(&header)?;
        if let Some(ref mut rw) = reject_writer {
            rw.finish(&header)?;
        }

        // Finalize and write metrics
        self.finalize_metrics(&mut umi_metrics, &unmatched_umi)?;

        // Log summary
        let rejected = missing_umis + wrong_length + mismatched;
        let kept = total_records - rejected;
        info!("Read {total_records}; kept {kept} and rejected {rejected}");
        info!("Total templates processed: {total_templates}");

        if missing_umis > 0 || wrong_length > 0 {
            warn!("###################################################################");
            if missing_umis > 0 {
                warn!("# {missing_umis} were missing UMI attributes in the BAM file!");
            }
            if wrong_length > 0 {
                warn!("# {wrong_length} had unexpected UMIs of differing lengths in the BAM file!");
            }
            warn!("###################################################################");
        }

        // Check minimum correction ratio
        if let Some(min) = self.min_corrected {
            #[allow(clippy::cast_precision_loss)]
            let ratio_kept = kept as f64 / total_records as f64;
            if ratio_kept < min {
                bail!(
                    "Final ratio of reads kept / total was {ratio_kept:.2} (user specified minimum was {min:.2}). \
                    This could indicate a mismatch between library preparation and the provided UMI file."
                );
            }
        }

        Ok(total_records)
    }
}

/// Counts mismatches between two sequences, stopping early if max is exceeded.
///
/// Efficiently counts mismatches between two byte sequences, stopping as soon as
/// `max_mismatches + 1` is reached. This early stopping is crucial for performance
/// when processing many sequences.
///
/// # Arguments
///
/// * `a` - First sequence
/// * `b` - Second sequence
/// * `max_mismatches` - Stop counting after this many mismatches
///
/// # Returns
///
/// The number of mismatches, capped at `max_mismatches + 1`.
///
/// # Example
///
/// ```
/// # use fgumi_lib::correct_umis::count_mismatches_with_max;
/// assert_eq!(count_mismatches_with_max(b"AAAAAA", b"AAAAAA", 10), 0);
/// assert_eq!(count_mismatches_with_max(b"AAAAAA", b"AAAAAT", 10), 1);
/// assert_eq!(count_mismatches_with_max(b"AAAAAA", b"CCCCCC", 2), 3);
/// ```
#[must_use]
pub fn count_mismatches_with_max(a: &[u8], b: &[u8], max_mismatches: usize) -> usize {
    let mut mismatches = 0;
    let min_len = a.len().min(b.len());

    for i in 0..min_len {
        if a[i] != b[i] {
            mismatches += 1;
            if mismatches > max_mismatches {
                return mismatches;
            }
        }
    }

    // Account for length differences
    mismatches += a.len().abs_diff(b.len());
    mismatches
}

/// A set of UMIs with pre-computed bit encodings for fast comparison.
///
/// Stores both the original byte sequences and their `BitEnc` representations
/// for efficient Hamming distance computation.
#[derive(Clone)]
pub struct EncodedUmiSet {
    /// Original byte sequences for fallback comparison
    bytes: Vec<Vec<u8>>,
    /// Bit-encoded sequences for fast comparison (None if encoding failed)
    encoded: Vec<Option<BitEnc>>,
    /// Original strings for output
    strings: Vec<String>,
}

impl EncodedUmiSet {
    /// Create a new `EncodedUmiSet` from string sequences.
    ///
    /// All sequences are converted to uppercase for case-insensitive matching.
    /// Pre-encodes all sequences that contain only ACGT bases.
    /// Sequences with other characters will use fallback byte comparison.
    pub fn new(sequences: &[String]) -> Self {
        // Store uppercase bytes for comparison
        let bytes: Vec<Vec<u8>> =
            sequences.iter().map(|s| s.bytes().map(|b| b.to_ascii_uppercase()).collect()).collect();
        let encoded: Vec<Option<BitEnc>> = bytes.iter().map(|b| BitEnc::from_bytes(b)).collect();
        // Store uppercase strings for output
        let strings: Vec<String> = sequences.iter().map(|s| s.to_uppercase()).collect();

        Self { bytes, encoded, strings }
    }

    /// Get the number of UMIs in the set.
    #[inline]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Check if the set is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }
}

/// Find the best matching UMI using bit-encoded comparison.
///
/// Uses fast XOR + popcount comparison when both observed and expected UMIs
/// can be bit-encoded. Falls back to byte comparison otherwise.
///
/// This function tracks the best match index during iteration and only clones
/// the matching string once at the end, avoiding repeated string allocations.
fn find_best_match_encoded(
    observed: &[u8],
    umi_set: &EncodedUmiSet,
    max_mismatches: usize,
    min_distance_diff: usize,
) -> UmiMatch {
    let mut best_index: Option<usize> = None;
    let mut best_mismatches = usize::MAX;
    let mut second_best_mismatches = usize::MAX;

    // Try to encode the observed UMI
    let observed_encoded = BitEnc::from_bytes(observed);

    match observed_encoded {
        Some(obs_enc) => {
            // Fast path: use bit-encoded comparison
            for (i, fixed_enc) in umi_set.encoded.iter().enumerate() {
                let mismatches = if let Some(enc) = fixed_enc {
                    // Both can be compared with BitEnc
                    enc.hamming_distance(&obs_enc) as usize
                } else {
                    // Fixed UMI couldn't be encoded, fall back to bytes
                    count_mismatches_with_max(observed, &umi_set.bytes[i], second_best_mismatches)
                };

                if mismatches < best_mismatches {
                    second_best_mismatches = best_mismatches;
                    best_mismatches = mismatches;
                    best_index = Some(i);
                } else if mismatches < second_best_mismatches {
                    second_best_mismatches = mismatches;
                }
            }
        }
        None => {
            // Slow path: observed UMI can't be encoded, use byte comparison
            for (i, fixed_umi) in umi_set.bytes.iter().enumerate() {
                let mismatches =
                    count_mismatches_with_max(observed, fixed_umi, second_best_mismatches);

                if mismatches < best_mismatches {
                    second_best_mismatches = best_mismatches;
                    best_mismatches = mismatches;
                    best_index = Some(i);
                } else if mismatches < second_best_mismatches {
                    second_best_mismatches = mismatches;
                }
            }
        }
    }

    // Build result with single clone at end
    let Some(idx) = best_index else {
        return UmiMatch { matched: false, umi: String::new(), mismatches: usize::MAX };
    };

    let matched = if best_mismatches <= max_mismatches {
        let distance_to_second = second_best_mismatches.saturating_sub(best_mismatches);
        distance_to_second >= min_distance_diff
    } else {
        false
    };

    UmiMatch { matched, umi: umi_set.strings[idx].clone(), mismatches: best_mismatches }
}

/// Finds pairs of UMIs within a specified edit distance.
///
/// Identifies all pairs of UMIs that are within `distance` edits of each other,
/// which helps detect potentially ambiguous UMI sets.
///
/// # Arguments
///
/// * `umis` - Slice of UMI sequences to check
/// * `distance` - Maximum edit distance threshold
///
/// # Returns
///
/// A vector of tuples `(umi1, umi2, actual_distance)` for pairs within the threshold.
///
/// # Example
///
/// ```
/// # use fgumi_lib::correct_umis::find_umi_pairs_within_distance;
/// let umis = vec!["AAAA".to_string(), "AAAT".to_string()];
/// let pairs = find_umi_pairs_within_distance(&umis, 2);
/// assert_eq!(pairs.len(), 1);
/// ```
#[must_use]
pub fn find_umi_pairs_within_distance(
    umis: &[String],
    distance: usize,
) -> Vec<(String, String, usize)> {
    let mut pairs = Vec::new();
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            let d = count_mismatches_with_max(umis[i].as_bytes(), umis[j].as_bytes(), distance + 1);
            if d <= distance {
                pairs.push((umis[i].clone(), umis[j].clone(), d));
            }
        }
    }
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::{io::Write as IoWrite, path::Path};
    use tempfile::{NamedTempFile, TempDir};

    /// Helper struct for managing temporary test file paths.
    /// Keeps `TempDir` alive for the lifetime of the struct.
    struct TestPaths {
        dir: TempDir,
        pub output: PathBuf,
        pub rejects: PathBuf,
        pub metrics: PathBuf,
    }

    impl TestPaths {
        fn new() -> Result<Self> {
            let dir = TempDir::new()?;
            Ok(Self {
                output: dir.path().join("output.bam"),
                rejects: dir.path().join("rejects.bam"),
                metrics: dir.path().join("metrics.txt"),
                dir,
            })
        }

        /// Create a numbered output file path (e.g., "output1.bam", "output2.bam")
        fn output_n(&self, n: usize) -> PathBuf {
            self.dir.path().join(format!("output{n}.bam"))
        }
    }

    const FIXED_UMIS: &[&str] = &["AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT"];

    /// Converts the fixed UMI strings to byte vectors for testing.
    ///
    /// # Returns
    ///
    /// `EncodedUmiSet` containing the fixed UMI sequences
    fn get_encoded_umi_set() -> EncodedUmiSet {
        let strings: Vec<String> = FIXED_UMIS.iter().map(|s| (*s).to_string()).collect();
        EncodedUmiSet::new(&strings)
    }

    /// Helper to create a `CorrectUmis` with specified input/output paths
    #[allow(dead_code)]
    fn create_corrector_with_paths(input: PathBuf, output: PathBuf) -> CorrectUmis {
        CorrectUmis {
            io: BamIoOptions { input, output },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    #[test]
    fn test_find_best_match_perfect() {
        let umi_set = get_encoded_umi_set();

        let hit1 = find_best_match_encoded(b"AAAAAA", &umi_set, 2, 2);
        assert!(hit1.matched);
        assert_eq!(hit1.mismatches, 0);
        assert_eq!(hit1.umi, "AAAAAA");

        let hit2 = find_best_match_encoded(b"CCCCCC", &umi_set, 2, 2);
        assert!(hit2.matched);
        assert_eq!(hit2.mismatches, 0);
        assert_eq!(hit2.umi, "CCCCCC");
    }

    #[test]
    fn test_find_best_match_with_mismatches() {
        let umi_set = get_encoded_umi_set();

        let m1 = find_best_match_encoded(b"AAAAAA", &umi_set, 2, 2);
        assert!(m1.matched);
        assert_eq!(m1.umi, "AAAAAA");
        assert_eq!(m1.mismatches, 0);

        let m2 = find_best_match_encoded(b"AAAAAT", &umi_set, 2, 2);
        assert!(m2.matched);
        assert_eq!(m2.umi, "AAAAAA");
        assert_eq!(m2.mismatches, 1);

        let m3 = find_best_match_encoded(b"AAAACT", &umi_set, 2, 2);
        assert!(m3.matched);
        assert_eq!(m3.umi, "AAAAAA");
        assert_eq!(m3.mismatches, 2);

        let m4 = find_best_match_encoded(b"AAAGCT", &umi_set, 2, 2);
        assert!(!m4.matched);
        assert_eq!(m4.mismatches, 3);
    }

    #[test]
    fn test_no_match_when_umis_too_similar() {
        let umi_set = EncodedUmiSet::new(&["AAAG".to_string(), "AAAT".to_string()]);

        let m = find_best_match_encoded(b"AAAG", &umi_set, 2, 2);
        assert!(!m.matched);
    }

    #[test]
    fn test_match_with_many_mismatches() {
        let umi_set = get_encoded_umi_set();

        let m1 = find_best_match_encoded(b"AAAAAA", &umi_set, 3, 2);
        assert!(m1.matched);
        assert_eq!(m1.umi, "AAAAAA");
        assert_eq!(m1.mismatches, 0);

        let m2 = find_best_match_encoded(b"AAACGT", &umi_set, 3, 2);
        assert!(m2.matched);
        assert_eq!(m2.umi, "AAAAAA");
        assert_eq!(m2.mismatches, 3);

        let m3 = find_best_match_encoded(b"AAACCC", &umi_set, 3, 2);
        assert!(!m3.matched);
        assert_eq!(m3.mismatches, 3);

        let m4 = find_best_match_encoded(b"AAACCT", &umi_set, 3, 2);
        assert!(!m4.matched);
        assert_eq!(m4.mismatches, 3);
    }

    #[test]
    fn test_find_umi_pairs_within_distance_none() {
        let pairs = find_umi_pairs_within_distance(
            &["AAAA".to_string(), "TTTT".to_string(), "CCCC".to_string(), "GGGG".to_string()],
            2,
        );
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_find_umi_pairs_within_distance_some() {
        let umis = vec![
            "ACACAC".to_string(),
            "CTCTCT".to_string(),
            "GAGAGA".to_string(),
            "TGTGTG".to_string(),
            "ACAGAC".to_string(),
            "AGAGAG".to_string(),
        ];
        let pairs = find_umi_pairs_within_distance(&umis, 2);
        assert_eq!(pairs.len(), 2);
        assert!(pairs.contains(&("ACACAC".to_string(), "ACAGAC".to_string(), 1)));
        assert!(pairs.contains(&("ACAGAC".to_string(), "AGAGAG".to_string(), 2)));
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement_str("AAAAAA"), "TTTTTT");
        assert_eq!(reverse_complement_str("AAAAGA"), "TCTTTT");
        assert_eq!(reverse_complement_str("ACGT"), "ACGT");
        assert_eq!(reverse_complement_str("GGGGGG"), "CCCCCC");
    }

    #[test]
    fn test_count_mismatches_with_max() {
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"AAAAAA", 10), 0);
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"AAAAAT", 10), 1);
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"AAAATT", 10), 2);
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"CCCCCC", 10), 6);

        // Test early stopping
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"CCCCCC", 2), 3);
    }
    /// Creates a test BAM file with specified records and UMI tags.
    ///
    /// Generates a minimal BAM file with a single reference sequence ("chr1")
    /// and the provided records. Each record is assigned default values for
    /// sequence, quality scores, and mapping information.
    ///
    /// # Arguments
    ///
    /// * `records` - Vector of tuples containing (`read_name`, `optional_umi_string`)
    ///
    /// # Returns
    ///
    /// A `NamedTempFile` handle to the created BAM file. The file is automatically
    /// deleted when the handle is dropped.
    ///
    /// # Errors
    ///
    /// Returns an error if BAM file creation or writing fails.
    fn create_test_bam(records: Vec<(&str, Option<&str>)>) -> Result<NamedTempFile> {
        use fgumi_lib::sam::builder::RecordBuilder;
        use noodles::sam::header::record::value::map::Map;

        let temp_file = NamedTempFile::new()?;
        let path = temp_file.path().to_path_buf();

        // Create a minimal BAM with header
        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZero::new(1000).unwrap(),
                ),
            )
            .build();

        let mut writer = noodles::bam::io::writer::Builder.build_from_path(&path)?;
        writer.write_header(&header)?;

        for (name, umi) in records {
            let mut builder = RecordBuilder::new()
                .name(name)
                .sequence("AAAAAAAAAA")
                .qualities(&[40u8; 10]) // 'I' is Phred+33 = 40
                .reference_sequence_id(0)
                .alignment_start(1)
                .mapping_quality(60);

            if let Some(umi_seq) = umi {
                builder = builder.tag("RX", umi_seq);
            }

            let record = builder.build();
            writer.write_alignment_record(&header, &record)?;
        }

        writer.finish(&header)?;
        Ok(temp_file)
    }

    /// Reads all record names from a BAM file.
    ///
    /// Opens the BAM file and extracts the name field from each record,
    /// useful for verifying which records were written to output files.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the BAM file to read
    ///
    /// # Returns
    ///
    /// Vector of read names as strings. Empty strings are used for records without names.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or read.
    fn read_bam_record_names(path: &Path) -> Result<Vec<String>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let _header = reader.read_header()?;

        let mut names = Vec::new();
        for result in reader.records() {
            let record = result?;
            let name = record.name().map(std::string::ToString::to_string).unwrap_or_default();
            names.push(name);
        }
        Ok(names)
    }

    #[test]
    fn test_rejects_reads_without_umis() -> Result<()> {
        let input_file = create_test_bam(vec![("q1", None)])?;
        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions {
                input: input_file.path().to_path_buf(),
                output: paths.output.clone(),
            },
            rejects_opts: RejectsOptions { rejects: Some(paths.rejects.clone()) },
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        assert_eq!(read_bam_record_names(&paths.output)?.len(), 0);
        assert_eq!(read_bam_record_names(&paths.rejects)?.len(), 1);

        Ok(())
    }

    #[test]
    fn test_validation_different_length_umis() {
        let temp_input = NamedTempFile::new().unwrap();
        let temp_output = NamedTempFile::new().unwrap();

        let corrector = CorrectUmis {
            io: BamIoOptions {
                input: temp_input.path().to_path_buf(),
                output: temp_output.path().to_path_buf(),
            },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "CCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert!(corrector.execute("test").is_err());
    }

    #[test]
    fn test_rejects_incorrect_length_umis() -> Result<()> {
        let input = create_test_bam(vec![("q1", Some("ACGT"))])?;
        let dir = TempDir::new()?;
        let output = dir.path().join("output.bam");
        let rejects = dir.path().join("rejects.bam");

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(rejects.clone()) },
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        assert_eq!(read_bam_record_names(&output)?.len(), 0);
        assert_eq!(read_bam_record_names(&rejects)?.len(), 1);

        Ok(())
    }

    #[test]
    fn test_end_to_end_single_umis() -> Result<()> {
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA")),
            ("q2", Some("AAAAAA")),
            ("q3", Some("AATAAA")),
            ("q4", Some("AATTAA")),
            ("q5", Some("AAATGC")),
            ("q6", Some("CCCCCC")),
            ("q7", Some("GGGGGG")),
            ("q8", Some("TTTTTT")),
            ("q9", Some("GGGTTT")),
            ("q10", Some("AAACCC")),
        ])?;

        let dir = TempDir::new()?;
        let output = dir.path().join("output.bam");
        let rejects = dir.path().join("rejects.bam");
        let metrics = dir.path().join("metrics.txt");

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(rejects.clone()) },
            metrics: Some(metrics.clone()),
            max_mismatches: 3,
            min_distance_diff: 2,
            umis: vec![
                "AAAAAA".to_string(),
                "CCCCCC".to_string(),
                "GGGGGG".to_string(),
                "TTTTTT".to_string(),
            ],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let corrected_names = read_bam_record_names(&output)?;
        assert_eq!(corrected_names.len(), 8);
        assert!(corrected_names.contains(&"q1".to_string()));
        assert!(corrected_names.contains(&"q8".to_string()));

        let rejected_names = read_bam_record_names(&rejects)?;
        assert_eq!(rejected_names.len(), 2);
        assert!(rejected_names.contains(&"q9".to_string()));
        assert!(rejected_names.contains(&"q10".to_string()));

        // Check metrics
        let metrics_data = UmiCorrectionMetrics::read_metrics(&metrics)?;
        let aaaaaa = metrics_data.iter().find(|m| m.umi == "AAAAAA").unwrap();
        assert_eq!(aaaaaa.total_matches, 5);
        assert_eq!(aaaaaa.perfect_matches, 2);
        assert_eq!(aaaaaa.one_mismatch_matches, 1);
        assert_eq!(aaaaaa.two_mismatch_matches, 1);
        assert_eq!(aaaaaa.other_matches, 1);

        Ok(())
    }

    #[test]
    fn test_end_to_end_duplex_umis() -> Result<()> {
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA-CCCCCC")),
            ("q2", Some("AAACAA-CCCCAC")),
            ("q3", Some("AAAAAA-ACTACT")),
            ("q4", Some("GGGGGG-TTTTTT")),
            ("q5", Some("GCGCGC-TTTTTT")),
        ])?;

        let dir = TempDir::new()?;
        let output = dir.path().join("output.bam");
        let rejects = dir.path().join("rejects.bam");

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(rejects.clone()) },
            metrics: None,
            max_mismatches: 3,
            min_distance_diff: 2,
            umis: vec![
                "AAAAAA".to_string(),
                "CCCCCC".to_string(),
                "GGGGGG".to_string(),
                "TTTTTT".to_string(),
            ],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let corrected_names = read_bam_record_names(&output)?;
        assert_eq!(corrected_names.len(), 3);
        assert!(corrected_names.contains(&"q1".to_string()));
        assert!(corrected_names.contains(&"q4".to_string()));

        let rejected_names = read_bam_record_names(&rejects)?;
        assert_eq!(rejected_names.len(), 2);
        assert!(rejected_names.contains(&"q3".to_string()));

        Ok(())
    }

    #[test]
    fn test_umi_loading_from_file() -> Result<()> {
        let input = create_test_bam(vec![("q1", Some("AAAAAA")), ("q2", Some("AATAAA"))])?;

        // Create UMI file
        let mut umi_file = NamedTempFile::new()?;
        writeln!(umi_file, "AAAAAA")?;
        writeln!(umi_file, "CCCCCC")?;
        writeln!(umi_file, "GGGGGG")?;
        writeln!(umi_file, "TTTTTT")?;
        umi_file.flush()?;

        let dir = TempDir::new()?;
        let output1 = dir.path().join("output1.bam");
        let metrics1 = dir.path().join("metrics1.txt");
        let output2 = dir.path().join("output2.bam");
        let metrics2 = dir.path().join("metrics2.txt");

        // Run with command line UMIs
        let corrector1 = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output1.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: Some(metrics1.clone()),
            max_mismatches: 3,
            min_distance_diff: 2,
            umis: vec![
                "AAAAAA".to_string(),
                "CCCCCC".to_string(),
                "GGGGGG".to_string(),
                "TTTTTT".to_string(),
            ],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        corrector1.execute("test")?;

        // Run with UMI file
        let corrector2 = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output2.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: Some(metrics2.clone()),
            max_mismatches: 3,
            min_distance_diff: 2,
            umis: vec![],
            umi_files: vec![umi_file.path().to_path_buf()],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        corrector2.execute("test")?;

        // Compare metrics
        let m1 = UmiCorrectionMetrics::read_metrics(&metrics1)?;
        let m2 = UmiCorrectionMetrics::read_metrics(&metrics2)?;

        assert_eq!(m1.len(), m2.len());
        for (metric1, metric2) in m1.iter().zip(m2.iter()) {
            assert_eq!(metric1.umi, metric2.umi);
            assert_eq!(metric1.total_matches, metric2.total_matches);
        }

        Ok(())
    }

    #[test]
    fn test_original_umi_storage() -> Result<()> {
        let input =
            create_test_bam(vec![("exact", Some("AAAAAA")), ("correctable", Some("AAAAGA"))])?;

        let dir = TempDir::new()?;
        let output = dir.path().join("output.bam");

        // Test with original UMI storage (default)
        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "TTTTTT".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        corrector.execute("test")?;

        // Read output and check OX tag
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&output)?;
        let header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let record_buf = RecordBuf::try_from_alignment_record(&header, &record)?;
            let name = record_buf.name().map(std::string::ToString::to_string).unwrap_or_default();

            let ox_tag = Tag::ORIGINAL_UMI_BARCODE_SEQUENCE;
            if name == "exact" {
                assert!(record_buf.data().get(&ox_tag).is_none());
            } else if name == "correctable" {
                let ox_value = record_buf.data().get(&ox_tag);
                assert!(ox_value.is_some());
            }
        }

        // Test with original UMI storage disabled
        let output2 = dir.path().join("output2.bam");
        let corrector2 = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output2.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "TTTTTT".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: true,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        corrector2.execute("test")?;

        let mut reader2 = noodles::bam::io::reader::Builder.build_from_path(&output2)?;
        let header2 = reader2.read_header()?;

        for result in reader2.records() {
            let record = result?;
            let record_buf = RecordBuf::try_from_alignment_record(&header2, &record)?;
            assert!(record_buf.data().get(&Tag::ORIGINAL_UMI_BARCODE_SEQUENCE).is_none());
        }

        Ok(())
    }

    #[test]
    fn test_revcomp_option() -> Result<()> {
        let input = create_test_bam(vec![
            ("exact", Some(&reverse_complement_str("AAAAAA"))),
            ("correctable", Some(&reverse_complement_str("AAAAGA"))),
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "TTTTTT".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: true,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        corrector.execute("test")?;

        // Read output and verify
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(&paths.output)?;
        let header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let record_buf = RecordBuf::try_from_alignment_record(&header, &record)?;
            let name = record_buf.name().map(std::string::ToString::to_string).unwrap_or_default();

            if name == "correctable" {
                let ox_tag = Tag::ORIGINAL_UMI_BARCODE_SEQUENCE;
                let ox_value = record_buf.data().get(&ox_tag);
                assert!(ox_value.is_some());
                // Original should be "TCTTTT" (revcomp of "AAAAGA")
                if let Some(sam::alignment::record_buf::data::field::Value::String(s)) = ox_value {
                    assert_eq!(s.to_string(), "TCTTTT");
                }
            }
        }

        Ok(())
    }

    #[test]
    fn test_exact_match_only_mode() -> Result<()> {
        // Test with max_mismatches=0 (exact match only)
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA")), // Perfect match
            ("q2", Some("AAAAAB")), // 1 mismatch - should reject
            ("q3", Some("CCCCCC")), // Perfect match
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(paths.rejects.clone()) },
            metrics: None,
            max_mismatches: 0, // Exact match only
            min_distance_diff: 1,
            umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let corrected_names = read_bam_record_names(&paths.output)?;
        assert_eq!(corrected_names.len(), 2);
        assert!(corrected_names.contains(&"q1".to_string()));
        assert!(corrected_names.contains(&"q3".to_string()));

        let rejected_names = read_bam_record_names(&paths.rejects)?;
        assert_eq!(rejected_names.len(), 1);
        assert!(rejected_names.contains(&"q2".to_string()));

        Ok(())
    }

    #[test]
    fn test_all_reads_already_correct() -> Result<()> {
        // Test scenario where all reads have perfect UMI matches
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA")),
            ("q2", Some("CCCCCC")),
            ("q3", Some("GGGGGG")),
            ("q4", Some("TTTTTT")),
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(paths.rejects.clone()) },
            metrics: Some(paths.metrics.clone()),
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec![
                "AAAAAA".to_string(),
                "CCCCCC".to_string(),
                "GGGGGG".to_string(),
                "TTTTTT".to_string(),
            ],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let corrected_names = read_bam_record_names(&paths.output)?;
        assert_eq!(corrected_names.len(), 4);

        let rejected_names = read_bam_record_names(&paths.rejects)?;
        assert_eq!(rejected_names.len(), 0);

        // All should be perfect matches (excluding the unmatched "NNNNNN" entry)
        let metrics_data = UmiCorrectionMetrics::read_metrics(&paths.metrics)?;
        for metric in &metrics_data {
            if !metric.umi.starts_with('N') {
                assert_eq!(metric.perfect_matches, 1);
                assert_eq!(metric.one_mismatch_matches, 0);
            }
        }

        Ok(())
    }

    #[test]
    fn test_single_known_umi() -> Result<()> {
        // Test with only one known UMI
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA")),
            ("q2", Some("AAAAAB")),
            ("q3", Some("AAAACC")),
            ("q4", Some("CCCCCC")), // No match
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(paths.rejects.clone()) },
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string()], // Only one UMI
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let corrected_names = read_bam_record_names(&paths.output)?;
        assert_eq!(corrected_names.len(), 3);

        let rejected_names = read_bam_record_names(&paths.rejects)?;
        assert_eq!(rejected_names.len(), 1);
        assert!(rejected_names.contains(&"q4".to_string()));

        Ok(())
    }

    #[test]
    fn test_umis_with_n_bases() -> Result<()> {
        // Test UMIs containing N (ambiguous base)
        let input = create_test_bam(vec![
            ("q1", Some("ANAAAA")), // Has N
            ("q2", Some("AANAAA")), // Has N
            ("q3", Some("AAAAAA")), // Perfect
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(paths.rejects.clone()) },
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        // N bases count as mismatches
        let corrected_names = read_bam_record_names(&paths.output)?;
        assert_eq!(corrected_names.len(), 3);
        // All should match since N is within 2 mismatches

        Ok(())
    }

    #[test]
    fn test_very_long_umis() -> Result<()> {
        // Test with longer UMIs (20 bases)
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAAAAAAAAAAAAAAAA")), // Perfect
            ("q2", Some("AAAAAAAAAAAAAAAAAAAB")), // 1 mismatch
            ("q3", Some("CCCCCCCCCCCCCCCCCCCC")), // Perfect
            ("q4", Some("CCCCCCCCCCCCCCCCCCCT")), // 1 mismatch
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions { rejects: Some(paths.rejects.clone()) },
            metrics: None,
            max_mismatches: 1,
            min_distance_diff: 5, // Large distance for long UMIs
            umis: vec!["AAAAAAAAAAAAAAAAAAAA".to_string(), "CCCCCCCCCCCCCCCCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let corrected_names = read_bam_record_names(&paths.output)?;
        assert_eq!(corrected_names.len(), 4);

        Ok(())
    }

    #[test]
    fn test_min_corrected_threshold() -> Result<()> {
        // Test the min_corrected parameter (fraction of reads that must pass overall)
        // First test: 75% threshold should pass when 3/4 reads are correctable
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA")), // Perfect match - passes
            ("q2", Some("AAAAAB")), // 1 mismatch - passes
            ("q3", Some("AAAAAC")), // 1 mismatch - passes
            ("q4", Some("TTTTTT")), // No match with min_distance_diff=2 - fails
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 1,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: Some(0.75), // Require 75% of reads to pass
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Should succeed - 3/4 = 75% pass
        corrector.execute("test")?;
        let corrected_names = read_bam_record_names(&paths.output)?;
        assert_eq!(corrected_names.len(), 3);

        // Now test that it fails when ratio is too low
        let input2 = create_test_bam(vec![
            ("q1", Some("AAAAAA")), // Perfect match - passes
            ("q2", Some("TTTTTT")), // No match - fails
            ("q3", Some("GGGGGG")), // No match - fails
            ("q4", Some("CCCCCC")), // No match - fails
        ])?;

        let output2 = paths.output_n(2);

        let corrector2 = CorrectUmis {
            io: BamIoOptions { input: input2.path().to_path_buf(), output: output2.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 1,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: Some(0.75), // Require 75% but only 25% will pass
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Should fail - only 1/4 = 25% pass, which is below 75% threshold
        let result = corrector2.execute("test");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Final ratio of reads kept"));

        Ok(())
    }

    #[test]
    fn test_case_sensitivity() -> Result<()> {
        // Ensure UMIs are handled case-insensitively (converted to uppercase)
        let input = create_test_bam(vec![
            ("q1", Some("aaaaaa")), // Lowercase
            ("q2", Some("AaAaAa")), // Mixed case
            ("q3", Some("AAAAAA")), // Uppercase
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 0, // Exact match
            min_distance_diff: 1,
            umis: vec!["AAAAAA".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        // All should match after case normalization
        let corrected_names = read_bam_record_names(&paths.output)?;
        assert_eq!(corrected_names.len(), 3);

        Ok(())
    }

    #[test]
    fn test_output_order_preserved_with_multiple_threads() -> Result<()> {
        // Create a test BAM with many records in a specific order
        let mut records: Vec<(&str, Option<&str>)> = Vec::new();

        // Create 100 records with names that encode their position
        // Use different UMIs to ensure some need correction and some don't
        let umis = ["AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT"];
        for i in 0..100 {
            let name = format!("read_{i:04}");
            // Leak the string to get a &'static str for the test
            let name_static: &'static str = Box::leak(name.into_boxed_str());
            let umi = umis[i % 4];
            records.push((name_static, Some(umi)));
        }

        let input = create_test_bam(records)?;
        let dir = TempDir::new()?;
        let output = dir.path().join("output.bam");

        // Run with 4 threads to test parallel processing order
        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec![
                "AAAAAA".to_string(),
                "CCCCCC".to_string(),
                "GGGGGG".to_string(),
                "TTTTTT".to_string(),
            ],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::new(4), // Multiple threads
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        // Read output and verify order matches input
        let output_names = read_bam_record_names(&output)?;
        assert_eq!(output_names.len(), 100);

        // Verify records are in sequential order
        for (i, name) in output_names.iter().enumerate() {
            let expected = format!("read_{i:04}");
            assert_eq!(
                name, &expected,
                "Output order mismatch at index {i}: expected {expected}, got {name}",
            );
        }

        Ok(())
    }

    #[test]
    fn test_output_order_preserved_with_8_threads() -> Result<()> {
        // Test with more threads and more records for better coverage
        let mut records: Vec<(&str, Option<&str>)> = Vec::new();

        let umis = ["AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT"];
        for i in 0..500 {
            let name = format!("rec_{i:05}");
            let name_static: &'static str = Box::leak(name.into_boxed_str());
            let umi = umis[i % 4];
            records.push((name_static, Some(umi)));
        }

        let input = create_test_bam(records)?;
        let dir = TempDir::new()?;
        let output = dir.path().join("output.bam");

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: vec![
                "AAAAAA".to_string(),
                "CCCCCC".to_string(),
                "GGGGGG".to_string(),
                "TTTTTT".to_string(),
            ],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::new(8), // 8 threads
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let output_names = read_bam_record_names(&output)?;
        assert_eq!(output_names.len(), 500);

        for (i, name) in output_names.iter().enumerate() {
            let expected = format!("rec_{i:05}");
            assert_eq!(
                name, &expected,
                "Output order mismatch at index {i}: expected {expected}, got {name}",
            );
        }

        Ok(())
    }

    // ============================================================================
    // Tests for template-level UMI correction functions
    // ============================================================================

    /// Helper to create a `RecordBuf` with a name and optional RX tag
    fn create_record_with_umi(name: &str, umi: Option<&str>) -> RecordBuf {
        use fgumi_lib::sam::builder::RecordBuilder;

        let mut builder = RecordBuilder::new().name(name).sequence("ACGT").qualities(&[40u8; 4]); // 'I' is Phred+33 = 40

        if let Some(umi_seq) = umi {
            builder = builder.tag("RX", umi_seq);
        }

        builder.build()
    }

    #[test]
    fn test_extract_and_validate_template_umi_empty() {
        let records: Vec<RecordBuf> = vec![];
        let umi_tag = Tag::new(b'R', b'X');
        let result = CorrectUmis::extract_and_validate_template_umi(&records, umi_tag);
        assert!(result.is_none());
    }

    #[test]
    fn test_extract_and_validate_template_umi_single_record_with_umi() {
        let records = vec![create_record_with_umi("read1", Some("AAAAAA"))];
        let umi_tag = Tag::new(b'R', b'X');
        let result = CorrectUmis::extract_and_validate_template_umi(&records, umi_tag);
        assert_eq!(result, Some("AAAAAA".to_string()));
    }

    #[test]
    fn test_extract_and_validate_template_umi_single_record_without_umi() {
        let records = vec![create_record_with_umi("read1", None)];
        let umi_tag = Tag::new(b'R', b'X');
        let result = CorrectUmis::extract_and_validate_template_umi(&records, umi_tag);
        assert!(result.is_none());
    }

    #[test]
    fn test_extract_and_validate_template_umi_paired_matching() {
        let records = vec![
            create_record_with_umi("read1", Some("AAAAAA")),
            create_record_with_umi("read1", Some("AAAAAA")),
        ];
        let umi_tag = Tag::new(b'R', b'X');
        let result = CorrectUmis::extract_and_validate_template_umi(&records, umi_tag);
        assert_eq!(result, Some("AAAAAA".to_string()));
    }

    #[test]
    #[should_panic(expected = "mismatched UMIs")]
    fn test_extract_and_validate_template_umi_paired_mismatched() {
        let records = vec![
            create_record_with_umi("read1", Some("AAAAAA")),
            create_record_with_umi("read1", Some("CCCCCC")),
        ];
        let umi_tag = Tag::new(b'R', b'X');
        CorrectUmis::extract_and_validate_template_umi(&records, umi_tag);
    }

    #[test]
    #[should_panic(expected = "inconsistent UMI presence")]
    fn test_extract_and_validate_template_umi_inconsistent_presence() {
        let records = vec![
            create_record_with_umi("read1", Some("AAAAAA")),
            create_record_with_umi("read1", None),
        ];
        let umi_tag = Tag::new(b'R', b'X');
        CorrectUmis::extract_and_validate_template_umi(&records, umi_tag);
    }

    #[test]
    fn test_extract_and_validate_template_umi_paired_no_umi() {
        let records =
            vec![create_record_with_umi("read1", None), create_record_with_umi("read1", None)];
        let umi_tag = Tag::new(b'R', b'X');
        let result = CorrectUmis::extract_and_validate_template_umi(&records, umi_tag);
        assert!(result.is_none());
    }

    #[test]
    fn test_compute_template_correction_perfect_match() {
        let umi_set = get_encoded_umi_set();
        let mut cache = None;

        let correction = CorrectUmis::compute_template_correction(
            "AAAAAA", 6, false, 2, 2, &umi_set, &mut cache,
        );

        assert!(correction.matched);
        assert_eq!(correction.corrected_umi, Some("AAAAAA".to_string()));
        assert_eq!(correction.original_umi, "AAAAAA");
        assert!(!correction.needs_correction);
        assert!(!correction.has_mismatches);
        assert_eq!(correction.rejection_reason, RejectionReason::None);
    }

    #[test]
    fn test_compute_template_correction_with_mismatch() {
        let umi_set = get_encoded_umi_set();
        let mut cache = None;

        let correction = CorrectUmis::compute_template_correction(
            "AAAAAT", 6, false, 2, 2, &umi_set, &mut cache,
        );

        assert!(correction.matched);
        assert_eq!(correction.corrected_umi, Some("AAAAAA".to_string()));
        assert_eq!(correction.original_umi, "AAAAAT");
        assert!(correction.needs_correction);
        assert!(correction.has_mismatches);
        assert_eq!(correction.rejection_reason, RejectionReason::None);
    }

    #[test]
    fn test_compute_template_correction_wrong_length() {
        let umi_set = get_encoded_umi_set();
        let mut cache = None;

        let correction =
            CorrectUmis::compute_template_correction("AAAA", 6, false, 2, 2, &umi_set, &mut cache);

        assert!(!correction.matched);
        assert!(correction.corrected_umi.is_none());
        assert_eq!(correction.rejection_reason, RejectionReason::WrongLength);
    }

    #[test]
    fn test_compute_template_correction_too_many_mismatches() {
        let umi_set = get_encoded_umi_set();
        let mut cache = None;

        let correction = CorrectUmis::compute_template_correction(
            "AAAGGG", 6, false, 2, 2, &umi_set, &mut cache,
        );

        assert!(!correction.matched);
        assert!(correction.corrected_umi.is_none());
        assert_eq!(correction.rejection_reason, RejectionReason::Mismatched);
    }

    #[test]
    fn test_compute_template_correction_with_revcomp() {
        let umi_set = get_encoded_umi_set();
        let mut cache = None;

        // TTTTTT revcomp -> AAAAAA
        let correction =
            CorrectUmis::compute_template_correction("TTTTTT", 6, true, 2, 2, &umi_set, &mut cache);

        assert!(correction.matched);
        assert_eq!(correction.corrected_umi, Some("AAAAAA".to_string()));
        assert!(correction.needs_correction); // revcomp triggers correction
        assert!(!correction.has_mismatches);
    }

    #[test]
    fn test_compute_template_correction_dual_umi() {
        // Create a UMI set with 6-base UMIs for dual-index testing
        let umi_set = get_encoded_umi_set();
        let mut cache = None;

        // Dual UMI: AAAAAA-CCCCCC (both should match)
        let correction = CorrectUmis::compute_template_correction(
            "AAAAAA-CCCCCC",
            6,
            false,
            2,
            2,
            &umi_set,
            &mut cache,
        );

        assert!(correction.matched);
        assert_eq!(correction.corrected_umi, Some("AAAAAA-CCCCCC".to_string()));
        assert!(!correction.needs_correction);
    }

    #[test]
    fn test_apply_correction_to_record_no_correction_needed() {
        let record = create_record_with_umi("read1", Some("AAAAAA"));
        let correction = TemplateCorrection {
            matched: true,
            corrected_umi: Some("AAAAAA".to_string()),
            original_umi: "AAAAAA".to_string(),
            needs_correction: false,
            has_mismatches: false,
            matches: vec![],
            rejection_reason: RejectionReason::None,
        };

        let result = CorrectUmis::apply_correction_to_record(
            record,
            &correction,
            Tag::new(b'R', b'X'),
            false,
        );

        // UMI should remain unchanged
        let rx_tag = Tag::new(b'R', b'X');
        let umi = result.data().get(&rx_tag);
        assert!(umi.is_some());
    }

    #[test]
    fn test_apply_correction_to_record_with_correction() {
        let record = create_record_with_umi("read1", Some("AAAAAT"));
        let correction = TemplateCorrection {
            matched: true,
            corrected_umi: Some("AAAAAA".to_string()),
            original_umi: "AAAAAT".to_string(),
            needs_correction: true,
            has_mismatches: true,
            matches: vec![],
            rejection_reason: RejectionReason::None,
        };

        let result = CorrectUmis::apply_correction_to_record(
            record,
            &correction,
            Tag::new(b'R', b'X'),
            false,
        );

        // Check RX tag was updated
        let rx_tag = Tag::new(b'R', b'X');
        if let Some(sam::alignment::record_buf::data::field::Value::String(s)) =
            result.data().get(&rx_tag)
        {
            assert_eq!(String::from_utf8_lossy(s), "AAAAAA");
        } else {
            panic!("RX tag not found or wrong type");
        }

        // Check OX tag was added
        let ox_tag = Tag::ORIGINAL_UMI_BARCODE_SEQUENCE;
        if let Some(sam::alignment::record_buf::data::field::Value::String(s)) =
            result.data().get(&ox_tag)
        {
            assert_eq!(String::from_utf8_lossy(s), "AAAAAT");
        } else {
            panic!("OX tag not found or wrong type");
        }
    }

    #[test]
    fn test_apply_correction_to_record_dont_store_original() {
        let record = create_record_with_umi("read1", Some("AAAAAT"));
        let correction = TemplateCorrection {
            matched: true,
            corrected_umi: Some("AAAAAA".to_string()),
            original_umi: "AAAAAT".to_string(),
            needs_correction: true,
            has_mismatches: true,
            matches: vec![],
            rejection_reason: RejectionReason::None,
        };

        // dont_store_original_umis = true
        let result = CorrectUmis::apply_correction_to_record(
            record,
            &correction,
            Tag::new(b'R', b'X'),
            true,
        );

        // Check RX tag was updated
        let rx_tag = Tag::new(b'R', b'X');
        if let Some(sam::alignment::record_buf::data::field::Value::String(s)) =
            result.data().get(&rx_tag)
        {
            assert_eq!(String::from_utf8_lossy(s), "AAAAAA");
        } else {
            panic!("RX tag not found or wrong type");
        }

        // OX tag should NOT be added
        let ox_tag = Tag::ORIGINAL_UMI_BARCODE_SEQUENCE;
        assert!(result.data().get(&ox_tag).is_none());
    }

    // ============================================================================
    // Tests for metrics output - unmatched UMI row and zero-count rows
    // ============================================================================

    #[test]
    fn test_metrics_includes_unmatched_umi_row() -> Result<()> {
        // Test that uncorrectable reads are tracked in an "NNNNNN" row
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA")), // Perfect match
            ("q2", Some("GGGTTT")), // No match (too far from all UMIs)
            ("q3", Some("CCCCCC")), // Perfect match
            ("q4", Some("ACTGAC")), // No match
            ("q5", None),           // Missing UMI - should also count as unmatched
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: Some(paths.metrics.clone()),
            max_mismatches: 1, // Strict matching
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        // Read metrics and verify unmatched row exists
        let metrics_data = UmiCorrectionMetrics::read_metrics(&paths.metrics)?;

        // Should have 3 rows: AAAAAA, CCCCCC, and NNNNNN (unmatched)
        assert_eq!(metrics_data.len(), 3, "Expected 3 UMI rows in metrics");

        // Find the unmatched UMI row
        let unmatched = metrics_data.iter().find(|m| m.umi == "NNNNNN");
        assert!(unmatched.is_some(), "Unmatched UMI row (NNNNNN) not found in metrics");

        let unmatched = unmatched.unwrap();
        // 3 uncorrectable reads: q2, q4, q5
        assert_eq!(
            unmatched.total_matches, 3,
            "Unmatched row should have 3 total_matches for uncorrectable reads"
        );
        // All mismatch categories should be 0 for unmatched row (matches fgbio)
        assert_eq!(unmatched.perfect_matches, 0);
        assert_eq!(unmatched.one_mismatch_matches, 0);
        assert_eq!(unmatched.two_mismatch_matches, 0);
        assert_eq!(unmatched.other_matches, 0);

        Ok(())
    }

    #[test]
    fn test_metrics_includes_all_umi_rows_even_with_zero_counts() -> Result<()> {
        // Test that all fixed UMIs appear in metrics even if they have zero matches
        let input = create_test_bam(vec![
            ("q1", Some("AAAAAA")), // Only AAAAAA matches
            ("q2", Some("AAAAAA")),
        ])?;

        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: Some(paths.metrics.clone()),
            max_mismatches: 1,
            min_distance_diff: 2,
            umis: vec![
                "AAAAAA".to_string(),
                "CCCCCC".to_string(), // No reads match this
                "GGGGGG".to_string(), // No reads match this
            ],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        // Read metrics
        let metrics_data = UmiCorrectionMetrics::read_metrics(&paths.metrics)?;

        // Should have 4 rows: AAAAAA, CCCCCC, GGGGGG, and NNNNNN (even though some have 0 matches)
        assert_eq!(metrics_data.len(), 4, "Expected 4 UMI rows in metrics (including zeros)");

        // Verify all UMIs are present
        assert!(metrics_data.iter().any(|m| m.umi == "AAAAAA"), "AAAAAA should be in metrics");
        assert!(
            metrics_data.iter().any(|m| m.umi == "CCCCCC"),
            "CCCCCC should be in metrics (even with 0 matches)"
        );
        assert!(
            metrics_data.iter().any(|m| m.umi == "GGGGGG"),
            "GGGGGG should be in metrics (even with 0 matches)"
        );
        assert!(
            metrics_data.iter().any(|m| m.umi == "NNNNNN"),
            "NNNNNN (unmatched) should be in metrics"
        );

        // Verify counts
        let cccccc = metrics_data.iter().find(|m| m.umi == "CCCCCC").unwrap();
        assert_eq!(cccccc.total_matches, 0, "CCCCCC should have 0 total_matches");

        let gggggg = metrics_data.iter().find(|m| m.umi == "GGGGGG").unwrap();
        assert_eq!(gggggg.total_matches, 0, "GGGGGG should have 0 total_matches");

        let aaaaaa = metrics_data.iter().find(|m| m.umi == "AAAAAA").unwrap();
        assert_eq!(aaaaaa.total_matches, 2, "AAAAAA should have 2 total_matches");

        Ok(())
    }

    #[test]
    fn test_metrics_unmatched_row_with_multithreaded() -> Result<()> {
        // Test that unmatched UMI row is correct with multi-threaded processing
        // This specifically tests the bug fix in execute_threads_mode
        let mut records: Vec<(&str, Option<&str>)> = Vec::new();

        // Create a mix of correctable and uncorrectable reads
        for i in 0..50 {
            let name = format!("read_{i:03}");
            let name_static: &'static str = Box::leak(name.into_boxed_str());
            let umi = match i % 5 {
                0 => Some("AAAAAA"), // Perfect match
                1 => Some("CCCCCC"), // Perfect match
                2 => Some("AAAAAB"), // 1 mismatch to AAAAAA
                3 => Some("GGGTTT"), // No match
                4 => Some("ACTGAC"), // No match
                _ => unreachable!(),
            };
            records.push((name_static, umi));
        }

        let input = create_test_bam(records)?;
        let paths = TestPaths::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions { input: input.path().to_path_buf(), output: paths.output.clone() },
            rejects_opts: RejectsOptions::default(),
            metrics: Some(paths.metrics.clone()),
            max_mismatches: 1,
            min_distance_diff: 2,
            umis: vec!["AAAAAA".to_string(), "CCCCCC".to_string()],
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading: ThreadingOptions::new(4), // Multi-threaded!
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        // Read metrics
        let metrics_data = UmiCorrectionMetrics::read_metrics(&paths.metrics)?;

        // Should have 3 rows: AAAAAA, CCCCCC, and NNNNNN
        assert_eq!(metrics_data.len(), 3, "Expected 3 UMI rows in metrics");

        // Verify unmatched row
        let unmatched = metrics_data.iter().find(|m| m.umi == "NNNNNN").unwrap();
        // 20 uncorrectable reads (i % 5 == 3 or 4, so 10 + 10 = 20)
        assert_eq!(unmatched.total_matches, 20, "Unmatched row should have 20 total_matches");

        // Verify matched UMIs
        let aaaaaa = metrics_data.iter().find(|m| m.umi == "AAAAAA").unwrap();
        // i % 5 == 0: 10 perfect, i % 5 == 2: 10 with 1 mismatch = 20 total
        assert_eq!(aaaaaa.total_matches, 20, "AAAAAA should have 20 total_matches");
        assert_eq!(aaaaaa.perfect_matches, 10);
        assert_eq!(aaaaaa.one_mismatch_matches, 10);

        let cccccc = metrics_data.iter().find(|m| m.umi == "CCCCCC").unwrap();
        assert_eq!(cccccc.total_matches, 10, "CCCCCC should have 10 total_matches");
        assert_eq!(cccccc.perfect_matches, 10);

        Ok(())
    }

    /// Parameterized test for all threading modes.
    ///
    /// Tests:
    /// - `None`: Single-threaded fast path, no pipeline
    /// - `Some(1)`: Pipeline with 1 thread
    /// - `Some(2)`: Pipeline with 2 threads
    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::pipeline_1(ThreadingOptions::new(1))]
    #[case::pipeline_2(ThreadingOptions::new(2))]
    fn test_threading_modes(#[case] threading: ThreadingOptions) -> Result<()> {
        let input = create_test_bam(vec![
            ("read1", Some("AAAAAA")),
            ("read2", Some("AAAAAG")), // 1 mismatch from AAAAAA
            ("read3", Some("CCCCCC")),
        ])?;
        let output = NamedTempFile::new()?;

        let corrector = CorrectUmis {
            io: BamIoOptions {
                input: input.path().to_path_buf(),
                output: output.path().to_path_buf(),
            },
            rejects_opts: RejectsOptions::default(),
            metrics: None,
            max_mismatches: 2,
            min_distance_diff: 2,
            umis: FIXED_UMIS.iter().map(|s| (*s).to_string()).collect(),
            umi_files: vec![],
            umi_tag: "RX".to_string(),
            dont_store_original_umis: false,
            cache_size: 100_000,
            min_corrected: None,
            revcomp: false,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        corrector.execute("test")?;

        let names = read_bam_record_names(output.path())?;
        assert_eq!(names.len(), 3, "Should have 3 records");

        Ok(())
    }

    // ========================================================================
    // Raw-byte helper for correct tests
    // ========================================================================

    /// Build a minimal raw BAM record with a UMI tag for testing.
    fn make_raw_bam_for_correct(name: &[u8], flag: u16, umi: &[u8]) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8;
        let seq_len = 4usize;
        let cigar_ops: &[u32] = if (flag & fgumi_lib::sort::bam_fields::flags::UNMAPPED) == 0 {
            &[(seq_len as u32) << 4]
        } else {
            &[]
        };
        let n_cigar_op = cigar_ops.len() as u16;
        let seq_bytes = seq_len.div_ceil(2);
        let total = 32 + l_read_name as usize + cigar_ops.len() * 4 + seq_bytes + seq_len;
        let mut buf = vec![0u8; total];

        buf[0..4].copy_from_slice(&0i32.to_le_bytes());
        buf[4..8].copy_from_slice(&100i32.to_le_bytes());
        buf[8] = l_read_name;
        buf[9] = 30; // mapq
        buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        buf[14..16].copy_from_slice(&flag.to_le_bytes());
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes());

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        let cigar_start = name_start + l_read_name as usize;
        for (i, &op) in cigar_ops.iter().enumerate() {
            let off = cigar_start + i * 4;
            buf[off..off + 4].copy_from_slice(&op.to_le_bytes());
        }

        // Append UMI tag: RX:Z:<umi>\0
        buf.extend_from_slice(b"RXZ");
        buf.extend_from_slice(umi);
        buf.push(0);

        buf
    }

    // ========================================================================
    // extract_and_validate_template_umi_raw tests
    // ========================================================================

    #[test]
    fn test_extract_and_validate_template_umi_raw_single_record() {
        let raw = make_raw_bam_for_correct(
            b"rea",
            fgumi_lib::sort::bam_fields::flags::PAIRED
                | fgumi_lib::sort::bam_fields::flags::FIRST_SEGMENT,
            b"AAAAAA",
        );
        let result = CorrectUmis::extract_and_validate_template_umi_raw(&[raw], [b'R', b'X']);
        assert_eq!(result, Some("AAAAAA".to_string()));
    }

    #[test]
    fn test_extract_and_validate_template_umi_raw_matching_pair() {
        let r1 = make_raw_bam_for_correct(
            b"rea",
            fgumi_lib::sort::bam_fields::flags::PAIRED
                | fgumi_lib::sort::bam_fields::flags::FIRST_SEGMENT,
            b"ACGTAC",
        );
        let r2 = make_raw_bam_for_correct(
            b"rea",
            fgumi_lib::sort::bam_fields::flags::PAIRED
                | fgumi_lib::sort::bam_fields::flags::LAST_SEGMENT,
            b"ACGTAC",
        );
        let result = CorrectUmis::extract_and_validate_template_umi_raw(&[r1, r2], [b'R', b'X']);
        assert_eq!(result, Some("ACGTAC".to_string()));
    }

    #[test]
    #[should_panic(expected = "mismatched UMIs")]
    fn test_extract_and_validate_template_umi_raw_mismatch_panics() {
        let r1 = make_raw_bam_for_correct(
            b"rea",
            fgumi_lib::sort::bam_fields::flags::PAIRED
                | fgumi_lib::sort::bam_fields::flags::FIRST_SEGMENT,
            b"AAAAAA",
        );
        let r2 = make_raw_bam_for_correct(
            b"rea",
            fgumi_lib::sort::bam_fields::flags::PAIRED
                | fgumi_lib::sort::bam_fields::flags::LAST_SEGMENT,
            b"CCCCCC",
        );
        let _ = CorrectUmis::extract_and_validate_template_umi_raw(&[r1, r2], [b'R', b'X']);
    }

    #[test]
    fn test_extract_and_validate_template_umi_raw_empty() {
        let result = CorrectUmis::extract_and_validate_template_umi_raw(&[], [b'R', b'X']);
        assert!(result.is_none());
    }

    #[test]
    fn test_extract_and_validate_template_umi_raw_no_umi_tag() {
        // Record with no UMI tag at all
        let mut raw = vec![0u8; 42]; // minimal record
        raw[8] = 4; // l_read_name
        raw[14..16].copy_from_slice(
            &(fgumi_lib::sort::bam_fields::flags::PAIRED
                | fgumi_lib::sort::bam_fields::flags::FIRST_SEGMENT)
                .to_le_bytes(),
        );
        raw[16..20].copy_from_slice(&4u32.to_le_bytes()); // l_seq = 4
        raw[32..36].copy_from_slice(b"rea\0");
        let result = CorrectUmis::extract_and_validate_template_umi_raw(&[raw], [b'R', b'X']);
        assert!(result.is_none());
    }

    // ========================================================================
    // apply_correction_to_raw tests
    // ========================================================================

    #[test]
    fn test_apply_correction_to_raw_corrects_umi() {
        use fgumi_lib::sort::bam_fields;

        let mut raw = make_raw_bam_for_correct(
            b"rea",
            bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT,
            b"AAAAAG", // original UMI with 1 mismatch
        );

        let correction = TemplateCorrection {
            matched: true,
            corrected_umi: Some("AAAAAA".to_string()),
            original_umi: "AAAAAG".to_string(),
            needs_correction: true,
            has_mismatches: true,
            matches: vec![],
            rejection_reason: RejectionReason::None,
        };

        CorrectUmis::apply_correction_to_raw(&mut raw, &correction, [b'R', b'X'], false);

        // Verify corrected UMI
        let umi = bam_fields::find_string_tag_in_record(&raw, b"RX");
        assert_eq!(umi, Some(b"AAAAAA".as_ref()));

        // Verify OX tag with original UMI
        let ox = bam_fields::find_string_tag_in_record(&raw, b"OX");
        assert_eq!(ox, Some(b"AAAAAG".as_ref()));
    }

    #[test]
    fn test_apply_correction_to_raw_no_correction_needed() {
        use fgumi_lib::sort::bam_fields;

        let mut raw = make_raw_bam_for_correct(
            b"rea",
            bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT,
            b"AAAAAA",
        );
        let orig_len = raw.len();

        let correction = TemplateCorrection {
            matched: true,
            corrected_umi: Some("AAAAAA".to_string()),
            original_umi: "AAAAAA".to_string(),
            needs_correction: false,
            has_mismatches: false,
            matches: vec![],
            rejection_reason: RejectionReason::None,
        };

        CorrectUmis::apply_correction_to_raw(&mut raw, &correction, [b'R', b'X'], false);

        // Record should be unchanged
        assert_eq!(raw.len(), orig_len);
        let umi = bam_fields::find_string_tag_in_record(&raw, b"RX");
        assert_eq!(umi, Some(b"AAAAAA".as_ref()));
    }

    #[test]
    fn test_apply_correction_to_raw_dont_store_original() {
        use fgumi_lib::sort::bam_fields;

        let mut raw = make_raw_bam_for_correct(
            b"rea",
            bam_fields::flags::PAIRED | bam_fields::flags::FIRST_SEGMENT,
            b"AAAAAG",
        );

        let correction = TemplateCorrection {
            matched: true,
            corrected_umi: Some("AAAAAA".to_string()),
            original_umi: "AAAAAG".to_string(),
            needs_correction: true,
            has_mismatches: true,
            matches: vec![],
            rejection_reason: RejectionReason::None,
        };

        CorrectUmis::apply_correction_to_raw(&mut raw, &correction, [b'R', b'X'], true);

        // Corrected UMI should be present
        let umi = bam_fields::find_string_tag_in_record(&raw, b"RX");
        assert_eq!(umi, Some(b"AAAAAA".as_ref()));

        // OX tag should NOT be present (dont_store_original_umis = true)
        let ox = bam_fields::find_string_tag_in_record(&raw, b"OX");
        assert!(ox.is_none());
    }
}

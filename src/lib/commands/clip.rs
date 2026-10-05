//! `Clip` command implementation.
//!
//! Clips reads in a BAM file to remove overlapping portions of read pairs.
//! This is useful for variant calling to avoid double-counting evidence from
//! overlapping portions of paired reads.

use crate::clipper::{ClippingMode, RawRecordClipper};
use crate::metrics::clip::{ClipCounts, ClippingMetricsCollection};
use crate::sam::SamTag;
use crate::template::{InsertSizeEnd, compute_insert_size_from_ends};
use crate::validation::validate_file_exists;
use anyhow::Result;
use clap::Parser;
use fgumi_raw_bam::{RawRecord, RawRecordView};
use std::path::PathBuf;

use super::command::Command;
use super::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};

/// Clips reads in a BAM file to remove overlaps
#[derive(Parser, Debug, Clone)]
#[command(
    name = "clip",
    about = "\x1b[38;5;173m[POST-CONSENSUS]\x1b[0m \x1b[36mClip overlapping reads in BAM files\x1b[0m",
    long_about = r#"
Clips reads from the same template. Ensures that at least N bases are clipped from any end of the read (i.e.
R1 5' end, R1 3' end, R2 5' end, and R2 3' end). Optionally clips reads from the same template to eliminate overlap
between the reads. This ensures that downstream processes, particularly variant calling, cannot double-count
evidence from the same template when both reads span a variant site in the same template.

Clipping overlapping reads is only performed on FR read pairs, and is implemented by clipping approximately half
the overlapping bases from each read. By default hard clipping is performed; the mode may be changed with
--clipping-mode.

Secondary alignments and supplemental alignments are not clipped, but are passed through into the output.

In order to correctly clip reads by template and update mate information, the input BAM must be either
queryname sorted or query grouped. If your input BAM is not in an appropriate order the sort can be
done in streaming fashion with, for example:

  fgumi sort -i in.bam --order queryname | fgumi clip -i /dev/stdin ...

The output is written in the same order as the input. To produce coordinate-sorted output, pipe the
result through a separate `fgumi sort`.

Any existing NM, UQ and MD tags are repaired, and mate-pair information is updated.

Three clipping modes are supported:
1. `soft` - soft-clip the bases and qualities.
2. `soft-with-mask` - soft-clip and mask the bases and qualities (make bases Ns and qualities the minimum).
3. `hard` - hard-clip the bases and qualities.

The --upgrade-clipping parameter will convert all existing clipping in the input to the given more stringent mode:
from `soft` to either `soft-with-mask` or `hard`, and `soft-with-mask` to `hard`. In all other cases, clipping remains
the same prior to applying any other clipping criteria.
"#
)]
#[allow(clippy::struct_excessive_bools)]
pub struct Clip {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Reference FASTA file (required for tag regeneration)
    #[arg(short = 'r', long = "reference", alias = "ref", required = true)]
    pub reference: PathBuf,

    /// Clipping mode: soft, soft-with-mask, or hard
    #[arg(short = 'c', long = "clipping-mode", default_value_t = ClippingMode::Hard)]
    pub clipping_mode: ClippingMode,

    /// Clip overlapping read pairs
    #[arg(long = "clip-overlapping-reads", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub clip_overlapping_reads: bool,

    /// Clip reads that extend past their mate's start position
    #[arg(
        long = "clip-bases-past-mate",
        alias = "clip-extending-past-mate",
        value_name = "true|false",
        default_value = "false",
        num_args = 0..=1,
        default_missing_value = "true",
        action = clap::ArgAction::Set,
        value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true,
    )]
    pub clip_extending_past_mate: bool,

    /// Minimum bases to clip from 5' end of R1
    #[arg(long = "read-one-five-prime", default_value = "0")]
    pub read_one_five_prime: usize,

    /// Minimum bases to clip from 3' end of R1
    #[arg(long = "read-one-three-prime", default_value = "0")]
    pub read_one_three_prime: usize,

    /// Minimum bases to clip from 5' end of R2
    #[arg(long = "read-two-five-prime", default_value = "0")]
    pub read_two_five_prime: usize,

    /// Minimum bases to clip from 3' end of R2
    #[arg(long = "read-two-three-prime", default_value = "0")]
    pub read_two_three_prime: usize,

    /// Upgrade existing clipping to the specified clipping mode
    #[arg(short = 'H', long = "upgrade-clipping", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub upgrade_clipping: bool,

    /// Automatically clip extended attributes that match read length (base modification tags
    /// MM/ML/am/bm are kept in step with the clipped read separately; listed tags that are not
    /// per-base data for the read, such as RG, are never clipped)
    #[arg(short = 'a', long = "auto-clip-attributes", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub auto_clip_attributes: bool,

    /// Output file for clipping metrics (produced on the declarative chain, the
    /// only execution path, regardless of `--threads`)
    #[arg(short = 'm', long = "metrics")]
    pub metrics: Option<PathBuf>,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline stats options
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

// ============================================================================
// Types for 7-step pipeline processing
// ============================================================================

/// Length of the record's leading hard clip (0 if its CIGAR does not start with `H`).
fn leading_hard_clip(record: &RawRecord) -> usize {
    record
        .cigar_ops_typed()
        .take_while(|op| op.kind() == fgumi_raw_bam::CigarKind::HardClip)
        .map(|op| op.len() as usize)
        .sum()
}

/// What [`ClipParams::clip_template`] did to one template.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct ClipTemplateOutcome {
    /// Overlap clipping removed bases from the pair.
    pub(crate) overlap_clipped: bool,
    /// Mate-extension clipping removed bases from the pair.
    pub(crate) extend_clipped: bool,
    /// Records whose modification tags (`MM`/`ML`/`MN`/`am`/`bm`) could not be kept in step
    /// with the clipped SEQ and were removed.
    pub(crate) modification_tags_removed: u64,
}

/// A record's state before clipping, kept to bring its modification tags in step afterwards.
struct PreClip {
    seq: Vec<u8>,
    reverse: bool,
    unmapped: bool,
    leading_hard_clip: usize,
}

/// Per-template clipping configuration, decoupled from `&Clip`.
///
/// The declarative chain's `process_fn` runs in a `move` closure that cannot borrow `&self`,
/// so the per-template clipping decision is driven from this small `Copy` value rather than
/// from `&Clip`. Capturing it lets the chain's worker closure share the exact same
/// per-template clipping logic (`clip_template`/`clip_pair`/`clip_fragment`) that the unit
/// tests exercise directly, instead of maintaining two copies that can silently drift apart.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ClipParams {
    /// Upgrade existing clipping to the configured mode before applying new clipping.
    upgrade_clipping: bool,
    /// Clip the overlapping bases of an FR read pair.
    clip_overlapping_reads: bool,
    /// Clip bases of a read that extend past its mate's start.
    clip_extending_past_mate: bool,
    /// Minimum bases to clip from the 5' end of R1.
    read_one_five_prime: usize,
    /// Minimum bases to clip from the 3' end of R1.
    read_one_three_prime: usize,
    /// Minimum bases to clip from the 5' end of R2.
    read_two_five_prime: usize,
    /// Minimum bases to clip from the 3' end of R2.
    read_two_three_prime: usize,
}

impl ClipParams {
    /// Builds the per-template clip configuration from the parsed `Clip` command.
    pub(crate) fn from_clip(clip: &Clip) -> Self {
        Self {
            upgrade_clipping: clip.upgrade_clipping,
            clip_overlapping_reads: clip.clip_overlapping_reads,
            clip_extending_past_mate: clip.clip_extending_past_mate,
            read_one_five_prime: clip.read_one_five_prime,
            read_one_three_prime: clip.read_one_three_prime,
            read_two_five_prime: clip.read_two_five_prime,
            read_two_three_prime: clip.read_two_three_prime,
        }
    }

    /// Clips one template's primary reads and repairs mate information.
    ///
    /// Clips the primary R1/R2 pair (or a lone primary fragment) — located by SAM flags via
    /// [`find_primary_pair_indices`], never positional `records[0]`/`records[1]` — so templates
    /// carrying secondary/supplementary reads (len > 2) are clipped too, matching fgbio
    /// `ClipBam`'s `Template.r1`/`r2` handling. For a pair it then fixes full mate-pair info
    /// (mate coords/strand, mate-unmapped flag, MQ/MC, TLEN) via [`set_mate_info_raw`] — matching
    /// fgbio's `SamPairUtil.setMateInfo` — and repairs mate info on any supplementary alignments.
    ///
    /// When `--upgrade-clipping` is set, existing clipping is first upgraded on *every* read of
    /// the template — including secondary/supplementary alignments — before the primary pair is
    /// clipped, matching fgbio `ClipBam`'s `template.allReads.foreach(upgradeAllClipping)`
    /// (`ClipBam.scala:123`). Doing this template-wide pre-pass here (rather than per-primary
    /// inside `clip_pair`/`clip_fragment`) is what lets supplementary reads' clipping be
    /// upgraded too, and keeps both threading paths in lockstep since both go through this method.
    ///
    /// This is the single shared per-template implementation the chain's worker closure
    /// runs (and the unit tests exercise directly). When `metrics` is `Some`, per-read
    /// base-clip counts are accumulated into it: the chain passes the calling worker's
    /// `PerThreadAccumulator` slot when `--metrics` is set and `None` otherwise, relying
    /// solely on the returned per-template flags for its atomic summary counters.
    ///
    /// The clipped bases' methylation calls are dropped from `MM`/`ML` (and `am`/`bm`), whether
    /// `soft-with-mask` masked them or hard clipping removed them, so the modification tags keep
    /// describing SEQ; after hard clipping `MN`, when present, is the new SEQ length. Tags that
    /// cannot be kept in step are removed and counted in the returned outcome.
    ///
    /// Returns, as a [`ClipTemplateOutcome`], whether overlap and/or mate-extension clipping
    /// removed any bases from the pair. A lone fragment reports neither, as do
    /// the non-clipping shapes fgbio also passes through untouched — a lone primary R2 (no R1) and
    /// an empty / all-secondary-and-supplementary template.
    ///
    /// # Errors
    ///
    /// Returns an error if primary-pair detection or a clipping operation fails.
    pub(crate) fn clip_template(
        &self,
        records: &mut [RawRecord],
        clipper: &RawRecordClipper,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<ClipTemplateOutcome> {
        // `soft-with-mask` clipping rewrites the clipped bases to N, hard clipping removes them,
        // and unmapping a reverse-mapped read reverse-complements SEQ, so snapshot every record
        // with modification tags and afterwards bring MM/ML in step with the clipped SEQ.
        let pre_clip: Vec<Option<PreClip>> = records
            .iter()
            .map(|record| {
                fgumi_consensus::filter::has_modification_tags(record).then(|| {
                    let view = RawRecordView::new(record);
                    PreClip {
                        seq: view.sequence_vec(),
                        reverse: view.is_reverse(),
                        unmapped: view.is_unmapped(),
                        leading_hard_clip: leading_hard_clip(record),
                    }
                })
            })
            .collect();
        let (overlap_clipped, extend_clipped) =
            self.clip_template_records(records, clipper, metrics)?;
        let mut modification_tags_removed = 0;
        for (record, pre_clip) in records.iter_mut().zip(&pre_clip) {
            let Some(pre_clip) = pre_clip else { continue };
            let view = RawRecordView::new(record);
            // Bases removed from the start of SEQ show up as added leading hard clip, unless
            // clipping then unmapped the read and cleared its CIGAR: that offset is unknown.
            let removed_start = if view.l_seq() as usize == pre_clip.seq.len() {
                Some(0)
            } else if view.is_unmapped() && !pre_clip.unmapped {
                None
            } else {
                leading_hard_clip(record).checked_sub(pre_clip.leading_hard_clip)
            };
            if fgumi_consensus::filter::trim_clipped_modifications_raw(
                record.as_mut_vec(),
                &pre_clip.seq,
                pre_clip.reverse,
                removed_start,
            ) {
                modification_tags_removed += 1;
            }
        }
        Ok(ClipTemplateOutcome { overlap_clipped, extend_clipped, modification_tags_removed })
    }

    /// Upgrades and clips one template's records; see [`Self::clip_template`].
    fn clip_template_records(
        &self,
        records: &mut [RawRecord],
        clipper: &RawRecordClipper,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<(bool, bool)> {
        // Upgrade existing clipping on *every* read of the template first — including
        // secondary/supplementary alignments — matching fgbio ClipBam (ClipBam.scala:123) before
        // clipping the primary pair. The per-read `clip_pair`/`clip_fragment` helpers do not run
        // this whole-read upgrade, so this pre-pass is its sole site for both threading paths.
        // (Those helpers still upgrade existing clipping at the specific end they clip, as fgbio's
        // `clip{5,3}PrimeEndOfRead` and `clipOverlappingReads` do, regardless of this flag.)
        if self.upgrade_clipping {
            for record in records.iter_mut() {
                clipper.upgrade_all_clipping_raw(record)?;
            }
        }

        match find_primary_pair_indices(records)? {
            (Some(i1), Some(i2)) => {
                let [r1, r2] =
                    records.get_disjoint_mut([i1, i2]).expect("distinct primary indices");
                let outcome = self.clip_pair(clipper, r1, r2, metrics)?;
                set_mate_info_raw(r1, r2);
                fix_supplemental_mate_info(records, i1, i2);
                Ok(outcome)
            }
            (Some(i1), None) => {
                self.clip_fragment(clipper, &mut records[i1], metrics)?;
                Ok((false, false))
            }
            // A lone primary R2 (second-of-pair with no first-of-pair mate) — or an empty /
            // all-secondary-and-supplementary template — is deliberately left unclipped and passed
            // through untouched. This mirrors fgbio `ClipBam`'s `case _ => ()` (`ClipBam.scala:133`):
            // fgbio clips only the `(r1, r2)` pair and the `(r1, None)` fragment cases, so a
            // second-of-pair primary orphaned from its R1 is a malformed template that fgbio does not
            // clip. Matching that exactly keeps fgumi's output byte-identical to fgbio's rather than
            // clipping the R2 with fragment/read-one thresholds it should not receive.
            (None, _) => Ok((false, false)),
        }
    }

    /// Clips a fragment (unpaired) read.
    ///
    /// Applies clipping operations to a single fragment read, including:
    /// 1. Applying fixed-position 5' and 3' clipping
    /// 2. Updating metrics if a metrics collector is provided
    ///
    /// Clipping upgrades (`--upgrade-clipping`) are *not* performed here: the caller runs a
    /// template-wide pre-pass that upgrades clipping on every read of the template (including
    /// secondary/supplementary alignments) before this method runs, matching fgbio `ClipBam`
    /// (`ClipBam.scala:123`).
    ///
    /// Fragment reads are treated as R1 for the purposes of fixed-position clipping.
    ///
    /// # Arguments
    ///
    /// * `clipper` - The record clipper instance
    /// * `record` - The fragment record to clip (mutable)
    /// * `metrics` - Optional metrics collector to update
    ///
    /// # Returns
    ///
    /// `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns an error if clipping operations fail.
    fn clip_fragment(
        &self,
        clipper: &RawRecordClipper,
        record: &mut RawRecord,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<()> {
        let prior_bases_clipped = clipped_bases_raw(record);

        // Note: clipping upgrades (`--upgrade-clipping`) are applied once per template over all
        // reads by the caller (matching fgbio ClipBam.scala:123), not here.

        // Apply fixed-position clipping
        let num_five_prime = if self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(record, self.read_one_five_prime)
        } else {
            0
        };

        let num_three_prime = if self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(record, self.read_one_three_prime)
        } else {
            0
        };

        // Update metrics
        if let Some(metrics) = metrics {
            metrics.fragment.update_raw(
                record,
                ClipCounts {
                    prior: prior_bases_clipped,
                    five_prime: num_five_prime,
                    three_prime: num_three_prime,
                    ..ClipCounts::default()
                },
            );
        }

        Ok(())
    }

    /// Clips a pair of reads with comprehensive clipping logic.
    ///
    /// Applies multiple types of clipping to a read pair:
    /// 1. Fixed-position 5' and 3' clipping for each read
    /// 2. Overlap clipping to remove duplicate coverage
    /// 3. Mate-extension clipping to remove reads extending past mate start
    /// 4. Updating metrics for both reads
    ///
    /// Clipping upgrades (`--upgrade-clipping`) are *not* performed here: the caller runs a
    /// template-wide pre-pass that upgrades clipping on every read of the template (including
    /// secondary/supplementary alignments) before this method runs, matching fgbio `ClipBam`
    /// (`ClipBam.scala:123`).
    ///
    /// The method intelligently determines which read is R1 vs R2 based on SAM flags
    /// and applies the appropriate fixed-position clipping thresholds.
    ///
    /// # Arguments
    ///
    /// * `clipper` - The record clipper instance
    /// * `r1` - First read of the pair (mutable)
    /// * `r2` - Second read of the pair (mutable)
    /// * `metrics` - Optional metrics collector to update
    ///
    /// # Returns
    ///
    /// A tuple of `(overlap_clipped, extend_clipped)` booleans indicating whether
    /// overlap or mate-extension clipping was performed.
    ///
    /// # Errors
    ///
    /// Returns an error if clipping operations fail.
    fn clip_pair(
        &self,
        clipper: &RawRecordClipper,
        r1: &mut RawRecord,
        r2: &mut RawRecord,
        metrics: Option<&mut ClippingMetricsCollection>,
    ) -> Result<(bool, bool)> {
        let prior_bases_clipped_r1 = clipped_bases_raw(r1);
        let prior_bases_clipped_r2 = clipped_bases_raw(r2);

        // Note: clipping upgrades (`--upgrade-clipping`) are applied once per template over all
        // reads by the caller (matching fgbio ClipBam.scala:123), not here.

        // Determine read types (raw flags: bit 6 = first segment, bit 7 = last segment)
        let (is_r1_first, is_r2_last) = (r1.is_first_segment(), r2.is_last_segment());

        // Apply fixed-position clipping for R1
        let num_r1_five_prime = if is_r1_first && self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r1, self.read_one_five_prime)
        } else if !is_r1_first && self.read_two_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r1, self.read_two_five_prime)
        } else {
            0
        };

        let num_r1_three_prime = if is_r1_first && self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r1, self.read_one_three_prime)
        } else if !is_r1_first && self.read_two_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r1, self.read_two_three_prime)
        } else {
            0
        };

        // Apply fixed-position clipping for R2
        let num_r2_five_prime = if is_r2_last && self.read_two_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r2, self.read_two_five_prime)
        } else if !is_r2_last && self.read_one_five_prime > 0 {
            clipper.clip_5_prime_end_of_read_raw(r2, self.read_one_five_prime)
        } else {
            0
        };

        let num_r2_three_prime = if is_r2_last && self.read_two_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r2, self.read_two_three_prime)
        } else if !is_r2_last && self.read_one_three_prime > 0 {
            clipper.clip_3_prime_end_of_read_raw(r2, self.read_one_three_prime)
        } else {
            0
        };

        // Clip overlapping reads
        let (num_overlapping_r1, num_overlapping_r2) = if self.clip_overlapping_reads {
            clipper.clip_overlapping_reads(r1, r2)
        } else {
            (0, 0)
        };

        // Clip reads extending past mate
        let (num_extending_r1, num_extending_r2) = if self.clip_extending_past_mate {
            clipper.clip_extending_past_mate_ends(r1, r2)
        } else {
            (0, 0)
        };

        // Update metrics
        if let Some(metrics) = metrics {
            let r1_counts = ClipCounts {
                prior: prior_bases_clipped_r1,
                five_prime: num_r1_five_prime,
                three_prime: num_r1_three_prime,
                overlapping: num_overlapping_r1,
                extending: num_extending_r1,
            };
            let r2_counts = ClipCounts {
                prior: prior_bases_clipped_r2,
                five_prime: num_r2_five_prime,
                three_prime: num_r2_three_prime,
                overlapping: num_overlapping_r2,
                extending: num_extending_r2,
            };

            // Determine which metric to update based on read flags
            if is_r1_first {
                metrics.read_one.update_raw(r1, r1_counts);
            } else {
                metrics.read_two.update_raw(r1, r1_counts);
            }

            if is_r2_last {
                metrics.read_two.update_raw(r2, r2_counts);
            } else {
                metrics.read_one.update_raw(r2, r2_counts);
            }
        }

        let overlap_clipped = num_overlapping_r1 > 0 || num_overlapping_r2 > 0;
        let extend_clipped = num_extending_r1 > 0 || num_extending_r2 > 0;

        Ok((overlap_clipped, extend_clipped))
    }
}

impl Command for Clip {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Reject two outputs resolving to one destination before any writer opens.
        let mut outputs: Vec<(&std::path::Path, &str)> =
            vec![(self.io.output.as_path(), "--output")];
        if let Some(path) = &self.metrics {
            outputs.push((path.as_path(), "--metrics"));
        }
        crate::commands::common::reject_output_collisions(&outputs)?;

        // Validate the input exists (stdin paths are exempt).
        self.io.validate()?;
        validate_file_exists(&self.reference, "Reference FASTA")?;

        // Validate clipping parameters (reader-free). At least one clipping option
        // must be requested. (`add_clip` re-checks this on the chain, but keeping it
        // here reports the error before any reader/writer opens.)
        if self.upgrade_clipping
            || self.clip_overlapping_reads
            || self.clip_extending_past_mate
            || self.read_one_five_prime > 0
            || self.read_one_three_prime > 0
            || self.read_two_five_prime > 0
            || self.read_two_three_prime > 0
        {
            // At least one clipping option is active
        } else {
            anyhow::bail!("At least one clipping option is required");
        }

        // The declarative chain is the only execution path. `execute` does the
        // reader-free pre-flight above (output collisions, including --metrics;
        // input + reference existence; clipping-option validation) and then always
        // dispatches to the chain, with or without `--threads` (absent `--threads`
        // runs the chain at a single worker). All user-facing diagnostics — the
        // `Clip` banner + Input/Output/mode lines, the `OperationTimer`, the
        // threading log lines, `require_query_grouped`, the summary counters, and
        // the `--metrics` TSV — are emitted inside `ChainBuilder::add_clip` and its
        // finalize hooks; running any of those here first would double-log and
        // pre-consume stdin. `--metrics` is produced on the chain too (#915).
        self.execute_chain(command_line)
    }
}

impl Clip {
    /// Runs the clip stage on the declarative chain builder
    /// (`ChainSpec::single_stage(Stage::Clip, ...)` → `build_for(spec)?.run()`).
    ///
    /// The chain is the only execution path: `execute` always dispatches here, with
    /// or without `--threads` (absent `--threads` runs the chain at a single
    /// worker). `add_clip` opens its own source, emits the timer/banner/threading
    /// log lines, loads the reference, validates the clipping options, and enforces
    /// `require_query_grouped`, so none of those run here — only the CRC-verify
    /// status line, which `add_clip` does not emit. `--metrics` is produced by the
    /// chain: `add_clip` builds a per-thread `ClippingMetricsCollection` accumulator
    /// when `self.metrics.is_some()` and reduces it in a success-only finalize hook.
    fn execute_chain(&self, command_line: &str) -> Result<()> {
        use crate::pipeline::chains::{
            ChainSpec, SingleStageContext, Stage, StageOptionsBag, build_for,
        };
        // add_clip emits the timer/banner/threading lines but NOT the CRC-verify
        // status line, so surface the effective CRC policy here.
        self.io.log_effective_check_crc();
        let stage_opts = StageOptionsBag { clip: Some(self.clone()), ..Default::default() };
        let ctx = SingleStageContext {
            io: &self.io,
            threading: &self.threading,
            compression: &self.compression,
            scheduler: &self.scheduler_opts,
            queue_memory: &self.queue_memory,
            command_line,
        };
        let spec = ChainSpec::single_stage(Stage::Clip, stage_opts, &ctx);
        build_for(spec)?.run()
    }
}

/// Snapshot of the mate-relevant fields of a read, taken before any mutation.
struct MateSnap {
    ref_id: i32,
    pos: i32,
    neg: bool,
    unmapped: bool,
    mapq: u8,
    cigar: String,
    aln_start: Option<usize>,
    aln_end: Option<usize>,
}

impl MateSnap {
    /// Snapshots only the fields `compute_insert_size_raw` reads, skipping the `cigar_to_string`
    /// allocation. `cigar` and `mapq` are left empty/zero, so the result must not be used to write
    /// a mate's MC or MQ tag.
    fn coords_of(rec: &RawRecord) -> Self {
        use fgumi_raw_bam::flags as rflags;
        Self {
            ref_id: rec.ref_id(),
            pos: rec.pos(),
            neg: rec.flags() & rflags::REVERSE != 0,
            unmapped: rec.flags() & rflags::UNMAPPED != 0,
            mapq: 0,
            cigar: String::new(),
            aln_start: rec.alignment_start_1based(),
            aln_end: rec.alignment_end_1based(),
        }
    }

    fn of(rec: &RawRecord) -> Self {
        use fgumi_raw_bam::flags as rflags;
        Self {
            ref_id: rec.ref_id(),
            pos: rec.pos(),
            neg: rec.flags() & rflags::REVERSE != 0,
            unmapped: rec.flags() & rflags::UNMAPPED != 0,
            mapq: rec.mapq(),
            cigar: rec.cigar_to_string(),
            aln_start: rec.alignment_start_1based(),
            aln_end: rec.alignment_end_1based(),
        }
    }
}

/// Sets or clears the `MATE_REVERSE` and `MATE_UNMAPPED` flags on a record.
fn set_mate_flags_raw(rec: &mut RawRecord, mate_neg: bool, mate_unmapped: bool) {
    use fgumi_raw_bam::flags as rflags;
    let mut f = rec.flags();
    if mate_neg {
        f |= rflags::MATE_REVERSE;
    } else {
        f &= !rflags::MATE_REVERSE;
    }
    if mate_unmapped {
        f |= rflags::MATE_UNMAPPED;
    } else {
        f &= !rflags::MATE_UNMAPPED;
    }
    rec.set_flags(f);
}

/// Writes the MQ (mate mapping quality) and MC (mate CIGAR) tags for a read.
fn set_mate_mq_mc_raw(rec: &mut RawRecord, mate_mapq: u8, mate_cigar: &str) {
    let mut editor = rec.tags_editor();
    editor.update_int(SamTag::MQ, i32::from(mate_mapq));
    editor.update_string(SamTag::MC, mate_cigar.as_bytes());
}

/// Removes the MQ and MC tags from a read (used when the mate is unmapped).
fn clear_mate_mq_mc_raw(rec: &mut RawRecord) {
    let mut editor = rec.tags_editor();
    editor.remove(SamTag::MQ);
    editor.remove(SamTag::MC);
}

/// Computes the TLEN (inferred insert size) for the first read of a pair, mirroring
/// htsjdk `SamPairUtil.computeInsertSize`. Returns 0 unless both reads are mapped to
/// the same reference; the second read's TLEN is the negation of this value.
///
/// Only the 5'-position extraction lives here — `MateSnap` caches the alignment ends, and a
/// read whose CIGAR yields none has no insert size. The arithmetic itself is
/// [`compute_insert_size_from_ends`], shared with `template`'s supplementary-TLEN fix-up so
/// the two commands cannot drift apart on ties, overflow, or the htsjdk adjustment.
fn compute_insert_size_raw(s1: &MateSnap, s2: &MateSnap) -> i32 {
    let five_prime = |s: &MateSnap| if s.neg { s.aln_end } else { s.aln_start };
    let (Some(p1), Some(p2)) = (five_prime(s1), five_prime(s2)) else {
        return 0;
    };
    let to_i64 = |p: usize| i64::try_from(p).unwrap_or(i64::MAX);
    compute_insert_size_from_ends(
        InsertSizeEnd::new(s1.ref_id, s1.unmapped, to_i64(p1)),
        InsertSizeEnd::new(s2.ref_id, s2.unmapped, to_i64(p2)),
    )
}

/// Sets full mate-pair information on a read pair, mirroring htsjdk
/// `SamPairUtil.setMateInfo(rec1, rec2, setMateCigar=true)`.
///
/// fgbio `ClipBam` calls `SamPairUtil.setMateInfo` on every pair after clipping so that
/// mate reference/position/strand, the mate-unmapped flag, the MQ/MC tags, and TLEN all
/// reflect the post-clip state — including reads that clipping unmapped (see CLIP-01). The
/// three branches (both mapped, both unmapped, one of each) match htsjdk exactly; an
/// unmapped read is relocated to its mapped mate's coordinate.
fn set_mate_info_raw(r1: &mut RawRecord, r2: &mut RawRecord) {
    // Snapshot both reads up front so mutating one never reads stale/updated fields of the other.
    let s1 = MateSnap::of(r1);
    let s2 = MateSnap::of(r2);

    if !s1.unmapped && !s2.unmapped {
        // Both mapped: copy each read's coordinates into the other's mate fields.
        r1.set_mate_ref_id(s2.ref_id);
        r1.set_mate_pos(s2.pos);
        set_mate_flags_raw(r1, s2.neg, false);
        set_mate_mq_mc_raw(r1, s2.mapq, &s2.cigar);

        r2.set_mate_ref_id(s1.ref_id);
        r2.set_mate_pos(s1.pos);
        set_mate_flags_raw(r2, s1.neg, false);
        set_mate_mq_mc_raw(r2, s1.mapq, &s1.cigar);

        let insert_size = compute_insert_size_raw(&s1, &s2);
        r1.set_template_length(insert_size);
        r2.set_template_length(-insert_size);
    } else if s1.unmapped && s2.unmapped {
        // Both unmapped: clear coordinates and mate coordinates, flag each mate unmapped.
        for (rec, mate_neg) in [(&mut *r1, s2.neg), (&mut *r2, s1.neg)] {
            rec.set_ref_id(-1);
            rec.set_pos(-1);
            rec.set_mate_ref_id(-1);
            rec.set_mate_pos(-1);
            set_mate_flags_raw(rec, mate_neg, true);
            clear_mate_mq_mc_raw(rec);
            rec.set_template_length(0);
            // Now unmapped (POS = -1): bin must be the SAM unmapped bin (4680).
            rec.recompute_bin();
        }
    } else {
        // Exactly one is unmapped: relocate it to the mapped mate's coordinate.
        let (mapped, unmapped, mapped_snap, unmapped_snap) = if s1.unmapped {
            (&mut *r2, &mut *r1, &s2, &s1)
        } else {
            (&mut *r1, &mut *r2, &s1, &s2)
        };

        // The unmapped read is placed at the mapped read's coordinate.
        unmapped.set_ref_id(mapped_snap.ref_id);
        unmapped.set_pos(mapped_snap.pos);
        unmapped.set_mate_ref_id(mapped_snap.ref_id);
        unmapped.set_mate_pos(mapped_snap.pos);
        set_mate_flags_raw(unmapped, mapped_snap.neg, false);
        set_mate_mq_mc_raw(unmapped, mapped_snap.mapq, &mapped_snap.cigar);
        unmapped.set_template_length(0);
        // POS moved from unmapped (-1) to the mate's coordinate; refresh the bin so
        // the placed read carries its position's bin (htsjdk recomputes on write).
        unmapped.recompute_bin();

        // The mapped read points its mate fields at the (now co-located) unmapped read.
        // Its mate-reverse flag reflects the unmapped read's *actual* strand (htsjdk
        // `SamPairUtil.java:267`): the clipper's own unmap clears REVERSE, but a read that
        // arrived unmapped-on-input may still carry it, and this runs on every pair.
        mapped.set_mate_ref_id(mapped_snap.ref_id);
        mapped.set_mate_pos(mapped_snap.pos);
        set_mate_flags_raw(mapped, unmapped_snap.neg, true);
        clear_mate_mq_mc_raw(mapped);
        mapped.set_template_length(0);
    }
}

/// Ports htsjdk 5.0.0 `SamPairUtil.setMateInformationOnSupplementalAlignment(supp, matePrimary,
/// setMateCigar=true)` **except for TLEN**, where it intentionally diverges.
///
/// htsjdk sets the supplementary's TLEN to the negation of the mate primary's. That assumes the
/// supplementary sits where its own primary sits, which is false by construction for split
/// alignments: the primary's TLEN describes coordinates this record does not occupy. Copying it
/// yields a non-zero TLEN across references, and the wrong sign and magnitude when the
/// supplementary lies beyond its mate. This computes TLEN from the supplementary's own alignment
/// against the mate primary instead, matching bwa-mem and minibwa. See issue #673 and
/// samtools/htsjdk#1795.
///
/// Everything else follows htsjdk: fgbio `ClipBam` calls this on every supplementary alignment
/// after clipping the primary pair, so a supplemental's mate fields point at its mate *primary*
/// read. It sets mate ref/pos/strand, the mate-unmapped flag, the mate CIGAR (MC, only when the
/// mate is mapped) and — as of htsjdk 5.0.0 — the mate mapping quality (MQ, unconditionally).
/// `mate` is snapshotted from the post-clip primary so the caller can fix several supplementals
/// without re-borrowing the primaries.
fn set_supplemental_mate_info_raw(supp: &mut RawRecord, mate: &MateSnap) {
    let tlen = compute_insert_size_raw(&MateSnap::coords_of(supp), mate);
    supp.set_mate_ref_id(mate.ref_id);
    supp.set_mate_pos(mate.pos);
    set_mate_flags_raw(supp, mate.neg, mate.unmapped);
    supp.set_template_length(tlen);
    let mut editor = supp.tags_editor();
    if mate.unmapped {
        editor.remove(SamTag::MC);
    } else {
        editor.update_string(SamTag::MC, mate.cigar.as_bytes());
    }
    // htsjdk 5.0.0 sets MQ unconditionally from the mate primary's mapping quality.
    editor.update_int(SamTag::MQ, i32::from(mate.mapq));
}

/// Finds the primary R1 and R2 record indices in a template, following fgbio `Template.r1`/`r2`:
/// the first non-secondary, non-supplementary read that is unpaired or first-of-pair, and the
/// first that is paired and second-of-pair, respectively. Either may be `None`.
///
/// Like fgbio's `Bams.Template` (`Bams.scala:161,167`), a template carrying more than one primary
/// (non-secondary, non-supplementary) R1 — or more than one primary R2 — is malformed and rejected
/// with an error rather than silently keeping the first and passing the extra through unclipped and
/// with unrepaired mate info. The error message mirrors fgbio verbatim.
fn find_primary_pair_indices(records: &[RawRecord]) -> Result<(Option<usize>, Option<usize>)> {
    let mut r1_idx = None;
    let mut r2_idx = None;
    for (i, rec) in records.iter().enumerate() {
        if rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        if !rec.is_paired() || rec.is_first_segment() {
            if r1_idx.is_some() {
                anyhow::bail!(
                    "Multiple non-secondary, non-supplemental R1s for {}",
                    String::from_utf8_lossy(rec.read_name()).trim_end_matches('\0')
                );
            }
            r1_idx = Some(i);
        } else if rec.is_last_segment() {
            if r2_idx.is_some() {
                anyhow::bail!(
                    "Multiple non-secondary, non-supplemental R2s for {}",
                    String::from_utf8_lossy(rec.read_name()).trim_end_matches('\0')
                );
            }
            r2_idx = Some(i);
        }
    }
    Ok((r1_idx, r2_idx))
}

/// Repairs mate information on the template's supplementary alignments after the primary pair
/// (at `r1_idx`/`r2_idx`) has been clipped and had its own mate info set. Mirrors fgbio
/// `ClipBam`: R1 supplementals point at the primary R2, R2 supplementals at the primary R1.
fn fix_supplemental_mate_info(records: &mut [RawRecord], r1_idx: usize, r2_idx: usize) {
    // Snapshot the post-clip primaries so the per-supplemental updates don't re-borrow them.
    let r1_snap = MateSnap::of(&records[r1_idx]);
    let r2_snap = MateSnap::of(&records[r2_idx]);

    for rec in records.iter_mut() {
        if !rec.is_supplementary() {
            continue;
        }
        // R1 supplementals (unpaired or first-of-pair) take R2 as their mate; R2 supplementals R1.
        if !rec.is_paired() || rec.is_first_segment() {
            set_supplemental_mate_info_raw(rec, &r2_snap);
        } else if rec.is_last_segment() {
            set_supplemental_mate_info_raw(rec, &r1_snap);
        }
    }
}

/// Returns the number of clipped bases (soft + hard) in a raw record's CIGAR.
fn clipped_bases_raw(record: &RawRecord) -> usize {
    record
        .cigar_ops_typed()
        .filter(|op| {
            matches!(
                op.kind(),
                fgumi_raw_bam::CigarKind::SoftClip | fgumi_raw_bam::CigarKind::HardClip
            )
        })
        .map(|op| op.len() as usize)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    // R2-CLIP-02: the primary R1/R2 are located by SAM flags, not by position, so templates
    // that carry secondary/supplementary alignments (len > 2) still resolve their primary pair.
    #[test]
    fn test_find_primary_pair_indices_ignores_secondary_and_supplementary() {
        use crate::sam::RecordBuilder;
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let enc = |b: &noodles::sam::alignment::RecordBuf| {
            encode_record_buf_to_raw(b, &header).expect("encode")
        };
        let mapped = |first: bool, secondary: bool, supplementary: bool, start: usize| {
            RecordBuilder::mapped_read()
                .name("q")
                .paired(true)
                .first_segment(first)
                .secondary(secondary)
                .supplementary(supplementary)
                .reference_sequence_id(0)
                .alignment_start(start)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build()
        };

        let r1 = mapped(true, false, false, 100);
        let r2 = mapped(false, false, false, 300);
        let supp_r1 = mapped(true, false, true, 700);
        let sec_r2 = mapped(false, true, false, 900);

        // Primaries at positions 0/1 with trailing secondary + supplementary reads.
        let recs = vec![enc(&r1), enc(&r2), enc(&supp_r1), enc(&sec_r2)];
        assert_eq!(find_primary_pair_indices(&recs).unwrap(), (Some(0), Some(1)));

        // Order-independent: a supplementary read first must not be mistaken for a primary.
        let recs2 = vec![enc(&supp_r1), enc(&r2), enc(&r1)];
        assert_eq!(find_primary_pair_indices(&recs2).unwrap(), (Some(2), Some(1)));

        // A lone fragment (unpaired primary) resolves R1 only.
        let frag = RecordBuilder::mapped_read()
            .name("f")
            .paired(false)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("50M")
            .sequence(&"A".repeat(50))
            .build();
        assert_eq!(find_primary_pair_indices(&[enc(&frag)]).unwrap(), (Some(0), None));
    }

    // A malformed template with two primary (non-secondary, non-supplementary) R1s — or two
    // primary R2s — is rejected loudly, matching fgbio `Bams.Template` (`Bams.scala:161,167`)
    // rather than silently keeping the first and passing the extra through unclipped.
    #[test]
    fn test_find_primary_pair_indices_rejects_duplicate_primaries() {
        use crate::sam::RecordBuilder;
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        let enc = |b: &noodles::sam::alignment::RecordBuf| {
            encode_record_buf_to_raw(b, &header).expect("encode")
        };
        let mapped = |first: bool, start: usize| {
            RecordBuilder::mapped_read()
                .name("dup")
                .paired(true)
                .first_segment(first)
                .reference_sequence_id(0)
                .alignment_start(start)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build()
        };

        // Two primary R1s (both first-of-pair, neither secondary/supplementary).
        let two_r1 =
            vec![enc(&mapped(true, 100)), enc(&mapped(false, 300)), enc(&mapped(true, 500))];
        let err = find_primary_pair_indices(&two_r1).unwrap_err().to_string();
        assert_eq!(err, "Multiple non-secondary, non-supplemental R1s for dup");

        // Two primary R2s (both last-of-pair).
        let two_r2 =
            vec![enc(&mapped(true, 100)), enc(&mapped(false, 300)), enc(&mapped(false, 500))];
        let err = find_primary_pair_indices(&two_r2).unwrap_err().to_string();
        assert_eq!(err, "Multiple non-secondary, non-supplemental R2s for dup");
    }

    /// Encodes a `RecordBuf` to a `RawRecord` with a shared single-contig header. Used by the
    /// raw mate-info tests below.
    fn encode_raw(rec: &noodles::sam::alignment::RecordBuf) -> RawRecord {
        use fgumi_raw_bam::encode_record_buf_to_raw;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        use std::num::NonZeroUsize;

        let ref_seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(100_000).expect("ref length must be nonzero"),
        );
        let header =
            noodles::sam::Header::builder().add_reference_sequence(b"chr1", ref_seq).build();
        encode_record_buf_to_raw(rec, &header).expect("encode")
    }

    /// Builds a mapped raw record from the common fields the mate-info tests need.
    fn raw_read(
        first: bool,
        supplementary: bool,
        reverse: bool,
        start: usize,
        mapq: u8,
    ) -> RawRecord {
        use crate::sam::RecordBuilder;
        encode_raw(
            &RecordBuilder::mapped_read()
                .name("q")
                .paired(true)
                .first_segment(first)
                .supplementary(supplementary)
                .reverse_complement(reverse)
                .reference_sequence_id(0)
                .alignment_start(start)
                .mapping_quality(mapq)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build(),
        )
    }

    // set_supplemental_mate_info_raw copies a mapped mate's coordinate/strand/MAPQ onto a
    // supplementary read, writes MC from the mate CIGAR, sets MQ, and computes TLEN from the
    // supplementary's own alignment against the mate.
    #[test]
    fn test_set_supplemental_mate_info_raw_mapped_mate() {
        use fgumi_raw_bam::flags as rflags;

        // A reverse-strand primary mate at 1-based 301 (0-based 300), MAPQ 40, CIGAR 50M.
        let mate = raw_read(false, false, true, 301, 40);
        let mut supp = raw_read(true, true, false, 700, 30);

        set_supplemental_mate_info_raw(&mut supp, &MateSnap::of(&mate));

        assert_eq!(supp.mate_ref_id(), 0);
        assert_eq!(supp.mate_pos(), 300);
        // The supplementary (700..749, forward, 5' = 700) is the rightmost segment; the mate
        // (301..350, reverse, 5' = 350) is leftmost. TLEN = 350 - 700 - 1.
        assert_eq!(supp.template_length(), -351);
        assert_ne!(supp.flags() & rflags::MATE_REVERSE, 0, "mate is reverse");
        assert_eq!(supp.flags() & rflags::MATE_UNMAPPED, 0, "mate is mapped");
        assert_eq!(supp.tags().find_mc(), Some(b"50M".as_slice()));
        assert_eq!(supp.tags().find_int(SamTag::MQ), Some(40));
    }

    // With an unmapped mate, set_supplemental_mate_info_raw flags the mate unmapped and drops MC,
    // but still writes MQ from the mate's mapping quality (htsjdk 5.0.0 sets MQ unconditionally).
    #[test]
    fn test_set_supplemental_mate_info_raw_unmapped_mate() {
        use crate::sam::RecordBuilder;
        use fgumi_raw_bam::flags as rflags;

        // Give the unmapped mate a distinct, non-zero MAPQ so the MQ assertion below proves the
        // value was copied from the mate rather than defaulting to 0.
        let mate = encode_raw(
            &RecordBuilder::mapped_read()
                .name("q")
                .paired(true)
                .first_segment(false)
                .unmapped(true)
                .reference_sequence_id(0)
                .alignment_start(301)
                .mapping_quality(37)
                .cigar("50M")
                .sequence(&"A".repeat(50))
                .build(),
        );
        // Seed an MC tag so we can confirm it is removed when the mate is unmapped.
        let mut supp = raw_read(true, true, false, 700, 30);
        supp.tags_editor().update_string(SamTag::MC, b"10M");
        assert!(supp.tags().contains(SamTag::MC), "MC present before");

        set_supplemental_mate_info_raw(&mut supp, &MateSnap::of(&mate));

        assert_eq!(supp.template_length(), 0, "TLEN is 0 when the mate is unmapped");
        assert_ne!(supp.flags() & rflags::MATE_UNMAPPED, 0, "mate is unmapped");
        assert!(!supp.tags().contains(SamTag::MC), "MC dropped when mate unmapped");
        // MQ is set unconditionally to the mate's mapping quality, even for an unmapped mate.
        assert_eq!(supp.tags().find_int(SamTag::MQ), Some(37), "MQ set from unmapped mate MAPQ");
    }

    // fix_supplemental_mate_info points R1 supplementals at the primary R2 and R2 supplementals
    // at the primary R1, inheriting each primary's coordinate/strand/MAPQ and computing TLEN from
    // the supplementary's own alignment.
    #[test]
    fn test_fix_supplemental_mate_info() {
        use fgumi_raw_bam::flags as rflags;

        let mut recs = vec![
            raw_read(true, false, false, 101, 60), // 0: primary R1 (forward, MAPQ 60)
            raw_read(false, false, true, 301, 40), // 1: primary R2 (reverse, MAPQ 40)
            raw_read(true, true, false, 701, 30),  // 2: supplementary R1
            raw_read(false, true, false, 901, 20), // 3: supplementary R2
        ];
        recs[0].set_template_length(200);
        recs[1].set_template_length(-200);

        fix_supplemental_mate_info(&mut recs, 0, 1);

        // Supp R1 (idx 2) takes primary R2 (idx 1) as its mate.
        assert_eq!(recs[2].mate_ref_id(), 0);
        assert_eq!(recs[2].mate_pos(), 300);
        assert_ne!(recs[2].flags() & rflags::MATE_REVERSE, 0, "primary R2 is reverse");
        assert_eq!(recs[2].tags().find_int(SamTag::MQ), Some(40));
        // Supp R1 (701..750, forward, 5' = 701) lies right of primary R2 (5' = 350): 350-701-1.
        assert_eq!(recs[2].template_length(), -352);

        // Supp R2 (idx 3) takes primary R1 (idx 0) as its mate.
        assert_eq!(recs[3].mate_ref_id(), 0);
        assert_eq!(recs[3].mate_pos(), 100);
        assert_eq!(recs[3].flags() & rflags::MATE_REVERSE, 0, "primary R1 is forward");
        assert_eq!(recs[3].tags().find_int(SamTag::MQ), Some(60));
        // Supp R2 (901..950, forward, 5' = 901) lies right of primary R1 (5' = 101): 101-901-1.
        assert_eq!(recs[3].template_length(), -801);
    }

    /// Builds a mapped record on an explicit reference for the supplementary-TLEN cases below.
    ///
    /// The builder validates `reference_sequence_id` against its single-contig test header, so the
    /// record is built on reference 0 and the id is stamped onto the encoded record afterwards.
    fn raw_read_on_ref(
        first: bool,
        supplementary: bool,
        reverse: bool,
        ref_id: i32,
        start: usize,
    ) -> RawRecord {
        let mut rec = raw_read(first, supplementary, reverse, start, 60);
        rec.set_ref_id(ref_id);
        rec
    }

    /// A supplementary's TLEN is computed from its own alignment against the mate primary, not
    /// copied from the mate primary's TLEN. See issue #673 and samtools/htsjdk#1795.
    ///
    /// The mate primary is always a reverse-strand read at 1-based 301 (50M, so 5' = 350).
    #[rstest]
    // Supplementary right of the mate: it is the rightmost segment, so TLEN is negative.
    #[case::beyond_mate(0, 700, -351)]
    // Supplementary left of the mate: it is the leftmost segment, so TLEN is positive.
    #[case::before_mate(0, 100, 251)]
    // Coincident 5' ends: htsjdk's convention gives the two ends differing signs via the +1/-1
    // adjustment, so the leftmost-by-tie-break gets +1.
    #[case::coincident_five_prime(0, 350, 1)]
    // Different reference: the information is unavailable, so TLEN is 0.
    #[case::cross_reference(1, 700, 0)]
    fn test_supplemental_tlen_is_computed_not_copied(
        #[case] supp_ref_id: i32,
        #[case] supp_start: usize,
        #[case] expected_tlen: i32,
    ) {
        let mate = raw_read_on_ref(false, false, true, 0, 301);
        let mut supp = raw_read_on_ref(true, true, false, supp_ref_id, supp_start);
        // Seed a value that the old copy-the-mate's-TLEN behaviour would have propagated, so a
        // regression cannot pass by coincidence.
        supp.set_template_length(-9999);

        set_supplemental_mate_info_raw(&mut supp, &MateSnap::of(&mate));

        assert_eq!(supp.template_length(), expected_tlen);
    }

    // The RecordBuf unsoftclipped_start/end helpers (used by the RecordBuf clip path) subtract or
    // add only *soft* clips, ignore hard clips, and return None for unmapped reads.
    #[test]
    fn test_unsoftclipped_recordbuf_helpers() {
        use crate::sam::RecordBuilder;
        use crate::sam::record_utils::{unsoftclipped_end, unsoftclipped_start};

        // 5H10S30M10S at 1-based 100: start = 100 - 10 (leading soft) = 90; hard clips ignored.
        // end = 100 + 30 (ref span) - 1 + 10 (trailing soft) = 139.
        let mapped = RecordBuilder::mapped_read()
            .name("q")
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("5H10S30M10S")
            .sequence(&"A".repeat(50))
            .build();
        assert_eq!(unsoftclipped_start(&mapped), Some(90));
        assert_eq!(unsoftclipped_end(&mapped), Some(139));

        let unmapped = RecordBuilder::mapped_read()
            .name("q")
            .unmapped(true)
            .reference_sequence_id(0)
            .alignment_start(100)
            .cigar("50M")
            .sequence(&"A".repeat(50))
            .build();
        assert_eq!(unsoftclipped_start(&unmapped), None);
        assert_eq!(unsoftclipped_end(&unmapped), None);
    }

    #[test]
    fn test_default_clip_parameters() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
        assert!(!clip.clip_overlapping_reads);
        assert!(!clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_with_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 5,
            read_one_three_prime: 3,
            read_two_five_prime: 7,
            read_two_three_prime: 2,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.read_one_five_prime, 5);
        assert_eq!(clip.read_one_three_prime, 3);
        assert_eq!(clip.read_two_five_prime, 7);
        assert_eq!(clip.read_two_three_prime, 2);
    }

    #[test]
    fn test_clip_with_overlapping_enabled() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_with_metrics_output() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: Some(PathBuf::from("metrics.txt")),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::SoftWithMask);
        assert!(clip.upgrade_clipping);
        assert_eq!(clip.metrics, Some(PathBuf::from("metrics.txt")));
    }

    #[test]
    fn test_clip_with_tag_regeneration() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.reference, PathBuf::from("reference.fa"));
        assert!(clip.auto_clip_attributes);
    }

    #[test]
    fn test_clip_all_modes_enabled() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 5,
            read_one_three_prime: 5,
            read_two_five_prime: 5,
            read_two_three_prime: 5,
            upgrade_clipping: true,
            auto_clip_attributes: true,
            metrics: Some(PathBuf::from("metrics.txt")),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // All options enabled
        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
        assert!(clip.upgrade_clipping);
        assert!(clip.auto_clip_attributes);
        assert!(clip.read_one_five_prime > 0);
    }

    #[test]
    fn test_clipping_mode_enum_values() {
        // Test that clipping_mode enum variants are set properly
        let soft = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(soft.clipping_mode, ClippingMode::Soft);
    }

    #[test]
    fn test_clip_asymmetric_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 10,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 15,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // R1 5' and R2 3' clipping only
        assert_eq!(clip.read_one_five_prime, 10);
        assert_eq!(clip.read_one_three_prime, 0);
        assert_eq!(clip.read_two_five_prime, 0);
        assert_eq!(clip.read_two_three_prime, 15);
    }

    #[test]
    fn test_clip_with_upgrade_all_clipping() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // upgrade_clipping should upgrade existing soft clips to hard clips
        assert!(clip.upgrade_clipping);
        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
    }

    #[test]
    fn test_clip_extending_past_mate_only() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Only clip_extending_past_mate is enabled
        assert!(!clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_overlapping_reads_only() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Only clip_overlapping_reads is enabled
        assert!(clip.clip_overlapping_reads);
        assert!(!clip.clip_extending_past_mate);
    }

    #[test]
    fn test_clip_modes_with_auto_clip_attributes() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // auto_clip_attributes should work with hard clipping
        assert!(clip.auto_clip_attributes);
        assert_eq!(clip.clipping_mode, ClippingMode::Hard);
    }

    #[test]
    fn test_clip_zero_bases_all_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // All fixed position clipping is zero (no fixed clipping)
        assert_eq!(clip.read_one_five_prime, 0);
        assert_eq!(clip.read_one_three_prime, 0);
        assert_eq!(clip.read_two_five_prime, 0);
        assert_eq!(clip.read_two_three_prime, 0);
    }

    #[test]
    fn test_clip_soft_with_mask_mode() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 5,
            read_one_three_prime: 5,
            read_two_five_prime: 5,
            read_two_three_prime: 5,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert_eq!(clip.clipping_mode, ClippingMode::SoftWithMask);
        assert!(clip.clip_overlapping_reads);
        assert_eq!(clip.read_one_five_prime, 5);
    }

    #[test]
    fn test_clip_large_fixed_positions() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 50,
            read_one_three_prime: 50,
            read_two_five_prime: 50,
            read_two_three_prime: 50,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Large fixed clipping values (e.g., for adapter trimming)
        assert_eq!(clip.read_one_five_prime, 50);
        assert_eq!(clip.read_one_three_prime, 50);
        assert_eq!(clip.read_two_five_prime, 50);
        assert_eq!(clip.read_two_three_prime, 50);
    }

    #[test]
    fn test_clip_combination_overlapping_and_fixed() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 10,
            read_one_three_prime: 10,
            read_two_five_prime: 10,
            read_two_three_prime: 10,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Both overlapping read clipping and fixed position clipping
        assert!(clip.clip_overlapping_reads);
        assert!(clip.read_one_five_prime > 0);
        assert!(clip.read_one_three_prime > 0);
    }

    #[test]
    fn test_clip_all_three_modes_comparison() {
        let soft = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        let soft_mask = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        let hard = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Verify all three modes are distinct
        assert_eq!(soft.clipping_mode, ClippingMode::Soft);
        assert_eq!(soft_mask.clipping_mode, ClippingMode::SoftWithMask);
        assert_eq!(hard.clipping_mode, ClippingMode::Hard);
    }

    #[test]
    fn test_clip_single_read_end_clipping() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 20,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        // Only R2 3' end clipping
        assert_eq!(clip.read_one_five_prime, 0);
        assert_eq!(clip.read_one_three_prime, 0);
        assert_eq!(clip.read_two_five_prime, 0);
        assert_eq!(clip.read_two_three_prime, 20);
    }

    #[test]
    fn test_clip_with_metrics_and_upgrade() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: Some(PathBuf::from("metrics.txt")),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert!(clip.clip_overlapping_reads);
        assert!(clip.upgrade_clipping);
        assert!(clip.metrics.is_some());
    }

    #[test]
    fn test_clip_both_extending_and_overlapping() {
        let clip = Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        assert!(clip.clip_overlapping_reads);
        assert!(clip.clip_extending_past_mate);
    }

    // Integration tests
    use crate::sam::{SamBuilder, Strand};
    use anyhow::Result;
    use tempfile::TempDir;

    fn create_test_reference(dir: &TempDir) -> PathBuf {
        let ref_path = dir.path().join("ref.fa");
        // Create a 200bp reference to accommodate read positions + padding
        let ref_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        std::fs::write(&ref_path, ref_content).expect("failed to write reference FASTA");
        // Also create index (200bp)
        let fai_content = "chr1\t200\t6\t200\t201\n";
        std::fs::write(dir.path().join("ref.fa.fai"), fai_content)
            .expect("failed to write FASTA index");
        ref_path
    }

    fn read_bam_records(path: &std::path::Path) -> Result<Vec<noodles::sam::alignment::RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let records: Vec<_> = reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?;
        Ok(records)
    }

    /// Renders a record's CIGAR as a string (e.g. `"5S15M"`) for assertions.
    fn cigar_string(record: &noodles::sam::alignment::RecordBuf) -> String {
        use noodles::sam::alignment::record::Cigar as _;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use std::fmt::Write as _;
        record.cigar().iter().filter_map(std::result::Result::ok).fold(
            String::new(),
            |mut acc, op| {
                let kind = match op.kind() {
                    Kind::Match => 'M',
                    Kind::Insertion => 'I',
                    Kind::Deletion => 'D',
                    Kind::Skip => 'N',
                    Kind::SoftClip => 'S',
                    Kind::HardClip => 'H',
                    Kind::Pad => 'P',
                    Kind::SequenceMatch => '=',
                    Kind::SequenceMismatch => 'X',
                };
                let _ = write!(acc, "{}{}", op.len(), kind);
                acc
            },
        )
    }

    /// Rewrites `builder`'s `@HD` line to advertise `SO:unsorted GO:query` (query grouped).
    /// `SamBuilder` writes a template's reads adjacently, so the emitted BAM is genuinely
    /// query grouped; the overlap-clip / clip-past-mate pipeline rejects input whose header
    /// does not advertise queryname sorting or query grouping (`clip` requires a template's
    /// reads to be adjacent).
    fn mark_query_grouped(builder: &mut SamBuilder) {
        builder.header = crate::sam::header_as_query_grouped(&builder.header);
    }

    /// CLIP3-02: fixed-position clipping must ensure *at least* N bases are clipped at
    /// the 5'/3' end *including any existing clipping*, matching fgbio `ClipBam`'s use of
    /// `clip5PrimeEndOfRead`/`clip3PrimeEndOfRead`. A read that already carries >= N
    /// bases of clipping at the requested end must not be clipped further.
    #[rstest]
    #[case::single_worker(ThreadingOptions::none())]
    #[case::multi_threaded(ThreadingOptions::new(2))]
    fn test_fixed_position_clip_counts_existing_clipping(
        #[case] threading: ThreadingOptions,
    ) -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // R1 forward with 5 bases already soft-clipped at its 5' (start) end.
        // R2 reverse with 4 bases already soft-clipped at its 5' (end) end.
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let _ = builder
            .add_pair()
            .name("read1")
            .contig(0)
            .start1(10)
            .cigar1("5S15M")
            .bases1("ACGTACGTACGTACGTACGT") // 20 bases
            .start2(120)
            .cigar2("16M4S")
            .bases2("ACGTACGTACGTACGTACGT") // 20 bases
            .build();
        mark_query_grouped(&mut builder);
        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            // Request fewer 5' bases than already exist on each read: fgbio clips 0 more there
            // (the existing-clipping check). Also request 5 bases at R1's 3' end, which has no
            // existing clipping, so clipping *is* applied there — this makes the command path
            // exercise real clipping rather than a no-op.
            read_one_five_prime: 3,
            read_one_three_prime: 5,
            read_two_five_prime: 2,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let records = read_bam_records(&output_path)?;
        assert_eq!(records.len(), 2);
        for record in &records {
            let cigar = cigar_string(record);
            if record.flags().is_first_segment() {
                assert_eq!(
                    cigar, "5S10M5S",
                    "R1: existing 5' clip counted (5S retained, no extra 5' clip), and the \
                     requested 5 bases newly clipped at the 3' end"
                );
            } else {
                assert_eq!(
                    cigar, "16M4S",
                    "R2 already had >= 2 bases 5'-clipped; expected no change"
                );
            }
        }
        Ok(())
    }

    /// CLIP3-01 (command level): overlap clipping must be strand-normalized end to end.
    /// A pair whose first-of-pair read is the reverse strand must yield the same
    /// per-strand output as the mirror pair whose first-of-pair read is the forward
    /// strand — the `clip` command must not clip the wrong (outer) ends. Exercised at
    /// both a single worker (no `--threads`) and a multi-worker pipeline.
    #[rstest]
    #[case::single_worker(0)]
    #[case::multi_threaded(4)]
    fn test_clip_overlap_negative_strand_first_matches_mirror(
        #[case] threads: usize,
    ) -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);

        // Runs `clip --clip-overlapping-reads` on an FR pair with the given strands and
        // returns (forward-read CIGAR, reverse-read CIGAR) from the output.
        let run = |tag: &str,
                   r1_minus: bool,
                   r1_start: usize,
                   r2_start: usize|
         -> Result<(String, String)> {
            let input_path = dir.path().join(format!("in_{tag}.bam"));
            let output_path = dir.path().join(format!("out_{tag}.bam"));
            let mut builder = SamBuilder::with_single_ref("chr1", 200);
            let _ = builder
                .add_pair()
                .name("read1")
                .contig(0)
                .start1(r1_start)
                .start2(r2_start)
                .strand1(if r1_minus { Strand::Minus } else { Strand::Plus })
                .strand2(if r1_minus { Strand::Plus } else { Strand::Minus })
                .build();
            mark_query_grouped(&mut builder);
            builder.write(&input_path)?;

            let threading = if threads == 0 {
                ThreadingOptions::none()
            } else {
                ThreadingOptions::new(threads)
            };
            let clip = Clip {
                io: BamIoOptions {
                    input: input_path,
                    output: output_path.clone(),
                    async_reader: false,
                    check_crc: false,
                    no_check_crc: false,
                },
                reference: ref_path.clone(),
                clipping_mode: ClippingMode::Soft,
                clip_overlapping_reads: true,
                clip_extending_past_mate: false,
                read_one_five_prime: 0,
                read_one_three_prime: 0,
                read_two_five_prime: 0,
                read_two_three_prime: 0,
                upgrade_clipping: false,
                auto_clip_attributes: false,
                metrics: None,
                threading,
                compression: CompressionOptions { compression_level: 1 },
                scheduler_opts: SchedulerOptions::default(),
                queue_memory: QueueMemoryOptions::default(),
            };
            clip.execute("test")?;

            let records = read_bam_records(&output_path)?;
            assert_eq!(records.len(), 2);
            let forward = records
                .iter()
                .find(|r| !r.flags().is_reverse_complemented())
                .map(cigar_string)
                .expect("a forward-strand read");
            let reverse = records
                .iter()
                .find(|r| r.flags().is_reverse_complemented())
                .map(cigar_string)
                .expect("a reverse-strand read");
            Ok((forward, reverse))
        };

        // Forward-first mirror: R1 forward on the left, R2 reverse on the right.
        let (fwd_forward, fwd_reverse) = run("fwdfirst", false, 10, 20)?;
        // Negative-strand first: R1 reverse on the right, R2 forward on the left.
        let (rev_forward, rev_reverse) = run("revfirst", true, 20, 10)?;

        // Pin the canonical (forward-first) output so the test fails if the command path
        // becomes a no-op: the 100bp reads overlap over [20, 110), so overlap clipping trims
        // 45bp from each read's 3' end where they meet (forward keeps [10, 65) -> 55M45S;
        // reverse keeps [65, 120) -> 45S55M). Both differ from the unclipped 100M input, so
        // a pipeline that silently skipped clipping would not satisfy these assertions.
        assert_eq!(fwd_forward, "55M45S", "forward-read must be 3'-clipped over the overlap");
        assert_eq!(fwd_reverse, "45S55M", "reverse-read must be 3'-clipped over the overlap");

        assert_eq!(
            fwd_forward, rev_forward,
            "forward-read CIGAR must not depend on R1/R2 strand order"
        );
        assert_eq!(
            fwd_reverse, rev_reverse,
            "reverse-read CIGAR must not depend on R1/R2 strand order"
        );
        Ok(())
    }

    /// CLIP3-03: `--upgrade-clipping` in `soft-with-mask` mode must mask existing
    /// soft-clipped bases to N with minimum quality while leaving the CIGAR intact,
    /// matching fgbio `ClipBam --upgrade-clipping --clipping-mode SoftWithMask`.
    #[test]
    fn test_upgrade_clipping_soft_with_mask_masks_existing_soft_bases() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // A fragment with 5 existing soft-clipped bases at the 5' end.
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        let _ = builder
            .add_frag()
            .name("maskme")
            .contig(0)
            .start(20)
            .cigar("5S15M")
            .bases("ACGTACGTACGTACGTACGT") // 20 bases
            .strand(Strand::Plus)
            .build();
        mark_query_grouped(&mut builder);
        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,
            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let records = read_bam_records(&output_path)?;
        assert_eq!(records.len(), 1);
        let record = &records[0];
        assert_eq!(cigar_string(record), "5S15M", "CIGAR unchanged");
        let bases = record.sequence().as_ref().to_vec();
        assert_eq!(&bases[..5], b"NNNNN", "leading soft bases masked to N");
        assert!(bases[5..].iter().all(|&b| b != b'N'), "aligned bases untouched");
        let quals = record.quality_scores().as_ref().to_vec();
        assert!(quals[..5].iter().all(|&q| q == 2), "leading soft quals masked to 2");
        Ok(())
    }

    #[test]
    fn test_clip_execute_basic_overlapping() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create overlapping read pair
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT") // 20 bases
            .contig(0)
            .start1(10) // R1 at pos 10
            .start2(20) // R2 at pos 20 - overlaps with R1
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert_eq!(output_records.len(), 2);

        Ok(())
    }

    #[test]
    fn test_clip_execute_soft_mode() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_soft_with_mask_mode() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::SoftWithMask,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_fixed_positions() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 3,
            read_one_three_prime: 2,
            read_two_five_prime: 2,
            read_two_three_prime: 3,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_extending_past_mate() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: true,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_upgrade_clipping() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: true,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_metrics() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let metrics_path = dir.path().join("metrics.txt");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: Some(metrics_path.clone()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());
        assert!(metrics_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_fragment() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create a fragment (unpaired) read
        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_frag()
            .name("frag1")
            .bases("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start(10)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 3,
            read_one_three_prime: 2,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert_eq!(output_records.len(), 1);

        Ok(())
    }

    #[test]
    fn test_clip_execute_with_auto_clip_attributes() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: true,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_all_clipping_options() -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let metrics_path = dir.path().join("metrics.txt");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();

        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: true,

            read_one_five_prime: 2,
            read_one_three_prime: 2,
            read_two_five_prime: 2,
            read_two_three_prime: 2,
            upgrade_clipping: true,
            auto_clip_attributes: true,
            metrics: Some(metrics_path.clone()),
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        assert!(output_path.exists());
        assert!(metrics_path.exists());

        Ok(())
    }

    #[test]
    fn test_clip_execute_no_clipping_option() {
        let dir = TempDir::new().expect("failed to create temp dir");
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(30)
            .build();

        builder.write(&input_path).expect("failed to write test BAM");

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path,
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false, // No clipping option
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };

        let result = clip.execute("test");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("At least one clipping option"));
    }

    /// Parameterized test for all threading modes (all run on the chain).
    ///
    /// Tests:
    /// - `None`: the chain at a single worker (no `--threads`)
    /// - `Some(1)`: pipeline with 1 thread
    /// - `Some(2)`: pipeline with 2 threads
    #[rstest]
    #[case::single_worker(ThreadingOptions::none())]
    #[case::pipeline_1(ThreadingOptions::new(1))]
    #[case::pipeline_2(ThreadingOptions::new(2))]
    fn test_threading_modes(#[case] threading: ThreadingOptions) -> Result<()> {
        let dir = TempDir::new()?;
        let ref_path = create_test_reference(&dir);
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut builder = SamBuilder::with_single_ref("chr1", 200);
        builder.set_queryname_sort_order(); // clip requires query-grouped input (CLIP3-05)
        let _ = builder
            .add_pair()
            .name("read1")
            .bases1("ACGTACGTACGTACGTACGT")
            .contig(0)
            .start1(10)
            .start2(20)
            .build();
        builder.write(&input_path)?;

        let clip = Clip {
            io: BamIoOptions {
                input: input_path,
                output: output_path.clone(),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: ref_path,
            clipping_mode: ClippingMode::Hard,
            clip_overlapping_reads: true,
            clip_extending_past_mate: false,

            read_one_five_prime: 0,
            read_one_three_prime: 0,
            read_two_five_prime: 0,
            read_two_three_prime: 0,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading,
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        };
        clip.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert_eq!(output_records.len(), 2, "Should have 2 records");

        Ok(())
    }

    /// Helper to create a `Clip` struct with specified clipping parameters and
    /// all other fields set to sensible defaults.
    fn make_clip(
        read_one_five_prime: usize,
        read_one_three_prime: usize,
        read_two_five_prime: usize,
        read_two_three_prime: usize,
    ) -> Clip {
        Clip {
            io: BamIoOptions {
                input: PathBuf::from("input.bam"),
                output: PathBuf::from("output.bam"),
                async_reader: false,
                check_crc: false,
                no_check_crc: false,
            },
            reference: PathBuf::from("reference.fa"),
            clipping_mode: ClippingMode::Soft,
            clip_overlapping_reads: false,
            clip_extending_past_mate: false,

            read_one_five_prime,
            read_one_three_prime,
            read_two_five_prime,
            read_two_three_prime,
            upgrade_clipping: false,
            auto_clip_attributes: false,
            metrics: None,
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    #[test]
    fn test_clip_fragment_with_metrics() {
        use crate::metrics::clip::ClippingMetricsCollection;

        let clip = make_clip(3, 2, 0, 0);
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        // Build a mapped fragment record with 20M CIGAR
        // SamBuilder pos is 0-based; alignment_start 100 → pos 99
        let mut record = fgumi_raw_bam::SamBuilder::new()
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .build();

        ClipParams::from_clip(&clip)
            .clip_fragment(&clipper, &mut record, Some(&mut metrics))
            .expect("clip_fragment should succeed");

        // Fragment metrics should be updated
        assert_eq!(metrics.fragment.reads, 1);
        assert_eq!(metrics.fragment.bases_clipped_five_prime, 3);
        assert_eq!(metrics.fragment.bases_clipped_three_prime, 2);
        // Remaining aligned bases: 20 - 3 - 2 = 15
        assert_eq!(metrics.fragment.bases, 15);
    }

    /// Hard clipping removes bases from SEQ, so the calls on them are dropped from `MM`/`ML`, the
    /// skips are recomputed over the remaining bases, and `MN` becomes the new SEQ length. `MM`
    /// indexes SEQ in original read orientation: on a reverse-mapped read the 5' bases are the
    /// end of SEQ as stored.
    #[rstest]
    #[case::forward(false, b"CACGTCAACG", b"C+m?,0,0,0;", &[10, 20, 30], b"CGTCAACG", b"C+m?,0,0;", &[20, 30])]
    #[case::reverse(true, b"CACGTCAACG", b"C+m?,0,0;", &[10, 20], b"CACGTCAA", b"C+m?,0;", &[20])]
    fn test_clip_template_hard_clipping_trims_modification_tags(
        #[case] reverse: bool,
        #[case] seq: &[u8],
        #[case] mm: &[u8],
        #[case] ml: &[u8],
        #[case] expected_seq: &[u8],
        #[case] expected_mm: &[u8],
        #[case] expected_ml: &[u8],
    ) {
        let clipper = RawRecordClipper::new(ClippingMode::Hard);
        let mut records = vec![
            fgumi_raw_bam::SamBuilder::new()
                .sequence(seq)
                .qualities(&[30; 10])
                .cigar_ops(&[10u32 << 4]) // 10M
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .flags(if reverse { fgumi_raw_bam::flags::REVERSE } else { 0 })
                .add_string_tag(SamTag::MM, mm)
                .add_array_u8(SamTag::ML, ml)
                .add_int_tag(SamTag::MN, 10)
                .build(),
        ];

        ClipParams::from_clip(&make_clip(2, 0, 0, 0))
            .clip_template(&mut records, &clipper, None)
            .expect("clip_template should succeed");

        let record = &records[0];
        assert_eq!(RawRecordView::new(record).sequence_vec(), expected_seq);
        let aux = fgumi_raw_bam::aux_data_slice(record);
        assert_eq!(fgumi_raw_bam::find_string_tag(aux, SamTag::MM), Some(expected_mm));
        assert_eq!(
            fgumi_raw_bam::find_array_tag(aux, SamTag::ML).map(|a| a.data.to_vec()).as_deref(),
            Some(expected_ml)
        );
        let expected_len = i64::try_from(expected_seq.len()).expect("short fixture");
        assert_eq!(fgumi_raw_bam::find_int_tag(aux, SamTag::MN), Some(expected_len));
    }

    /// One `clip_template` modification-tag case: the input record, the clipping, and the tags
    /// expected afterwards.
    struct ClipModCase {
        flags: u16,
        cigar: &'static [u32],
        seq: &'static [u8],
        mm: &'static str,
        ml: &'static [u8],
        five_prime: usize,
        three_prime: usize,
        upgrade_clipping: bool,
        auto_clip_attributes: bool,
        expected: ExpectedClipTags,
    }

    /// Expected SEQ, `MM`, `ML` and `MN` after clipping (`None` when the tags were removed).
    type ExpectedClipTags =
        (&'static [u8], Option<&'static str>, Option<&'static [u8]>, Option<i64>);

    const M10: &[u32] = &[10 << 4];

    /// Hard clipping keeps `MM`/`ML` in step with SEQ across the clipper's paths: an existing
    /// leading hard clip, upgraded soft clips, auto-clipped attributes, and reads that clipping
    /// unmaps. When clipping removes bases and then unmaps the read, the CIGAR no longer says
    /// where they came from, so the tags are removed rather than misplaced.
    #[rstest]
    #[case::existing_leading_hard_clip(ClipModCase {
        flags: 0, cigar: &[(3 << 4) | 5, 10 << 4], seq: b"CACGTCAACG", mm: "C+m?,0,0,0;",
        ml: &[10, 20, 30], five_prime: 5, three_prime: 0, upgrade_clipping: false,
        auto_clip_attributes: false,
        expected: (b"CGTCAACG", Some("C+m?,0,0;"), Some(&[20, 30][..]), Some(8)),
    })]
    #[case::upgraded_soft_clip(ClipModCase {
        flags: 0, cigar: &[(2 << 4) | 4, 8 << 4], seq: b"CACGTCAACG", mm: "C+m?,0,0,0;",
        ml: &[10, 20, 30], five_prime: 0, three_prime: 0, upgrade_clipping: true,
        auto_clip_attributes: false,
        expected: (b"CGTCAACG", Some("C+m?,0,0;"), Some(&[20, 30][..]), Some(8)),
    })]
    #[case::auto_clip_attributes_leave_mm_alone(ClipModCase {
        flags: 0, cigar: &[9 << 4], seq: b"CACGTCAAC", mm: "C+m?,0,0;", ml: &[10, 20],
        five_prime: 2, three_prime: 0, upgrade_clipping: false, auto_clip_attributes: true,
        expected: (b"CGTCAAC", Some("C+m?,0;"), Some(&[20][..]), Some(7)),
    })]
    #[case::reverse_read_unmapped(ClipModCase {
        flags: fgumi_raw_bam::flags::REVERSE, cigar: M10, seq: b"CACGTCAACG", mm: "C+m?,0,0;",
        ml: &[10, 20], five_prime: 10, three_prime: 0, upgrade_clipping: false,
        auto_clip_attributes: false,
        expected: (b"CGTTGACGTG", Some("C+m?,0,0;"), Some(&[10, 20][..]), Some(10)),
    })]
    #[case::hard_clipped_then_unmapped(ClipModCase {
        flags: 0, cigar: M10, seq: b"CACACACACA", mm: "C+m?,0,0,0,0,0;", ml: &[1, 2, 3, 4, 5],
        five_prime: 2, three_prime: 10, upgrade_clipping: false, auto_clip_attributes: false,
        expected: (b"CACACACA", None, None, None),
    })]
    #[case::leading_hard_clip_then_unmapped(ClipModCase {
        flags: 0, cigar: &[(5 << 4) | 5, 10 << 4], seq: b"CACGTCAACG", mm: "C+m?,0,0,0;",
        ml: &[10, 20, 30], five_prime: 15, three_prime: 0, upgrade_clipping: false,
        auto_clip_attributes: false,
        expected: (b"CACGTCAACG", Some("C+m?,0,0,0;"), Some(&[10, 20, 30][..]), Some(10)),
    })]
    fn test_clip_template_modification_tags_across_clipper_paths(#[case] case: ClipModCase) {
        let clipper =
            RawRecordClipper::with_auto_clip(ClippingMode::Hard, case.auto_clip_attributes);
        let mn = i32::try_from(case.seq.len()).expect("short fixture");
        let mut records = vec![
            fgumi_raw_bam::SamBuilder::new()
                .sequence(case.seq)
                .qualities(&vec![30; case.seq.len()])
                .cigar_ops(case.cigar)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .flags(case.flags)
                .add_string_tag(SamTag::MM, case.mm.as_bytes())
                .add_array_u8(SamTag::ML, case.ml)
                .add_int_tag(SamTag::MN, mn)
                .build(),
        ];
        let mut clip = make_clip(case.five_prime, case.three_prime, 0, 0);
        clip.upgrade_clipping = case.upgrade_clipping;

        let outcome = ClipParams::from_clip(&clip)
            .clip_template(&mut records, &clipper, None)
            .expect("clip_template should succeed");

        let (want_seq, want_mm, want_ml, want_mn) = case.expected;
        let record = &records[0];
        assert_eq!(RawRecordView::new(record).sequence_vec(), want_seq, "SEQ");
        let aux = fgumi_raw_bam::aux_data_slice(record);
        assert_eq!(
            fgumi_raw_bam::find_string_tag(aux, SamTag::MM),
            want_mm.map(str::as_bytes),
            "MM"
        );
        assert_eq!(
            fgumi_raw_bam::find_array_tag(aux, SamTag::ML).map(|a| a.data.to_vec()).as_deref(),
            want_ml,
            "ML"
        );
        assert_eq!(fgumi_raw_bam::find_int_tag(aux, SamTag::MN), want_mn, "MN");
        assert_eq!(outcome.modification_tags_removed, u64::from(want_mm.is_none()), "removed");
    }

    /// Clipping a read completely unmaps it, which reverse-complements a reverse-mapped read's
    /// SEQ into read orientation. `MM`/`ML` already index SEQ in that orientation and no base
    /// was masked, so the tags are kept as they are. `CACGTCAACG` reverse reads `CGTTGACGTG`.
    #[rstest]
    #[case::soft(ClippingMode::Soft)]
    #[case::soft_with_mask(ClippingMode::SoftWithMask)]
    fn test_clip_template_unmapping_keeps_modification_tags(#[case] mode: ClippingMode) {
        let clipper = RawRecordClipper::new(mode);
        let mut records = vec![
            fgumi_raw_bam::SamBuilder::new()
                .sequence(b"CACGTCAACG")
                .qualities(&[30; 10])
                .cigar_ops(&[10u32 << 4]) // 10M
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .flags(fgumi_raw_bam::flags::REVERSE)
                .add_string_tag(SamTag::MM, b"C+m?,0,0;")
                .add_array_u8(SamTag::ML, &[10, 20])
                .add_int_tag(SamTag::MN, 10)
                .build(),
        ];

        ClipParams::from_clip(&make_clip(10, 0, 0, 0))
            .clip_template(&mut records, &clipper, None)
            .expect("clip_template should succeed");

        let record = &records[0];
        assert!(RawRecordView::new(record).is_unmapped(), "clipped completely");
        let aux = fgumi_raw_bam::aux_data_slice(record);
        assert_eq!(fgumi_raw_bam::find_string_tag(aux, SamTag::MM), Some(&b"C+m?,0,0;"[..]));
        assert_eq!(
            fgumi_raw_bam::find_array_tag(aux, SamTag::ML).map(|a| a.data.to_vec()).as_deref(),
            Some(&[10u8, 20][..])
        );
    }

    /// `soft-with-mask` masks the clipped bases to N; the methylation calls on them are dropped
    /// from MM/ML so the tags keep describing SEQ. Plain soft clipping leaves SEQ and the tags
    /// alone.
    #[rstest]
    #[case::soft_with_mask(ClippingMode::SoftWithMask, b"NNCGTCAACG", b"C+m?,0,0;", &[20, 30])]
    #[case::soft(ClippingMode::Soft, b"CACGTCAACG", b"C+m?,0,0,0;", &[10, 20, 30])]
    fn test_clip_template_updates_modification_tags(
        #[case] mode: ClippingMode,
        #[case] expected_seq: &[u8],
        #[case] expected_mm: &[u8],
        #[case] expected_ml: &[u8],
    ) {
        let clipper = RawRecordClipper::new(mode);
        let mut records = vec![
            fgumi_raw_bam::SamBuilder::new()
                .sequence(b"CACGTCAACG")
                .qualities(&[30; 10])
                .cigar_ops(&[10u32 << 4]) // 10M
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .add_string_tag(SamTag::MM, b"C+m?,0,0,0;")
                .add_array_u8(SamTag::ML, &[10, 20, 30])
                .add_int_tag(SamTag::MN, 10)
                .build(),
        ];

        ClipParams::from_clip(&make_clip(2, 0, 0, 0))
            .clip_template(&mut records, &clipper, None)
            .expect("clip_template should succeed");

        let record = &records[0];
        assert_eq!(RawRecordView::new(record).sequence_vec(), expected_seq);
        let aux = fgumi_raw_bam::aux_data_slice(record);
        assert_eq!(fgumi_raw_bam::find_string_tag(aux, SamTag::MM), Some(expected_mm));
        assert_eq!(
            fgumi_raw_bam::find_array_tag(aux, SamTag::ML).map(|a| a.data.to_vec()).as_deref(),
            Some(expected_ml)
        );
    }

    #[test]
    fn test_clip_fragment_no_clipping_with_metrics() {
        use crate::metrics::clip::ClippingMetricsCollection;

        let clip = make_clip(0, 0, 0, 0);
        // clip_fragment only uses read_one_*; we just test the method directly.
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        let mut record = fgumi_raw_bam::SamBuilder::new()
            .sequence(b"ACGTACGTACGT")
            .qualities(&[30; 12])
            .cigar_ops(&[12u32 << 4]) // 12M
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .build();

        ClipParams::from_clip(&clip)
            .clip_fragment(&clipper, &mut record, Some(&mut metrics))
            .expect("clip_fragment should succeed");

        assert_eq!(metrics.fragment.reads, 1);
        assert_eq!(metrics.fragment.bases, 12);
        assert_eq!(metrics.fragment.bases_clipped_five_prime, 0);
        assert_eq!(metrics.fragment.bases_clipped_three_prime, 0);
    }

    #[test]
    fn test_clip_pair_with_metrics() {
        use crate::metrics::clip::ClippingMetricsCollection;
        use fgumi_raw_bam::flags as raw_flags;

        let clip = make_clip(2, 1, 1, 2);
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        // R1: first segment, forward strand, 20M
        // SamBuilder pos is 0-based; alignment_start 100 → pos 99
        let mut r1 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(199) // 0-based (mate_alignment_start 200)
            .template_length(120)
            .build();

        // R2: last segment, reverse strand, 20M (non-overlapping)
        let mut r2 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
            .ref_id(0)
            .pos(199) // 0-based (alignment_start 200)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(99) // 0-based (mate_alignment_start 100)
            .template_length(-120)
            .build();

        let (overlap, extend) = ClipParams::from_clip(&clip)
            .clip_pair(&clipper, &mut r1, &mut r2, Some(&mut metrics))
            .expect("clip_pair should succeed");

        assert!(!overlap, "no overlapping clipping expected");
        assert!(!extend, "no extending clipping expected");

        // R1 is first segment -> gets read_one clipping: 5'=2, 3'=1
        assert_eq!(metrics.read_one.reads, 1);
        assert_eq!(metrics.read_one.bases_clipped_five_prime, 2);
        assert_eq!(metrics.read_one.bases_clipped_three_prime, 1);

        // R2 is last segment -> gets read_two clipping: 5'=1, 3'=2
        assert_eq!(metrics.read_two.reads, 1);
        assert_eq!(metrics.read_two.bases_clipped_five_prime, 1);
        assert_eq!(metrics.read_two.bases_clipped_three_prime, 2);
    }

    #[test]
    fn test_clip_pair_with_metrics_swapped_flags() {
        use crate::metrics::clip::ClippingMetricsCollection;
        use fgumi_raw_bam::flags as raw_flags;

        // Test the !is_r1_first and !is_r2_last branches:
        // r1 is NOT first_segment, r2 is NOT last_segment
        let clip = make_clip(2, 1, 1, 2);
        let clipper = RawRecordClipper::new(ClippingMode::Soft);
        let mut metrics = ClippingMetricsCollection::new();

        // r1: last segment (not first)
        let mut r1 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT)
            .ref_id(0)
            .pos(99) // 0-based (alignment_start 100)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(199) // 0-based (mate_alignment_start 200)
            .template_length(120)
            .build();

        // r2: first segment (not last)
        let mut r2 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"pair1")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4]) // 20M
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::REVERSE)
            .ref_id(0)
            .pos(199) // 0-based (alignment_start 200)
            .mapq(60)
            .mate_ref_id(0)
            .mate_pos(99) // 0-based (mate_alignment_start 100)
            .template_length(-120)
            .build();

        let (overlap, extend) = ClipParams::from_clip(&clip)
            .clip_pair(&clipper, &mut r1, &mut r2, Some(&mut metrics))
            .expect("clip_pair should succeed");

        assert!(!overlap);
        assert!(!extend);

        // r1 is not first_segment -> goes to read_two metrics
        // r1 gets read_two clipping: 5'=1, 3'=2
        assert_eq!(metrics.read_two.reads, 1);
        assert_eq!(metrics.read_two.bases_clipped_five_prime, 1);
        assert_eq!(metrics.read_two.bases_clipped_three_prime, 2);

        // r2 is not last_segment -> goes to read_one metrics
        // r2 gets read_one clipping: 5'=2, 3'=1
        assert_eq!(metrics.read_one.reads, 1);
        assert_eq!(metrics.read_one.bases_clipped_five_prime, 2);
        assert_eq!(metrics.read_one.bases_clipped_three_prime, 1);
    }

    /// `set_mate_info_raw` mirrors htsjdk `SamPairUtil.setMateInfo` across all three
    /// branches (both mapped, one unmapped, both unmapped) so `fgumi clip` matches fgbio
    /// `ClipBam` when clipping unmaps a read (CLIP-01).
    #[test]
    fn test_set_mate_info_raw_all_branches() {
        use fgumi_raw_bam::flags as raw_flags;

        let has_flag = |rec: &RawRecord, f: u16| rec.flags() & f != 0;
        let has_tag = |rec: &RawRecord, tag: [u8; 2]| {
            fgumi_raw_bam::find_tag_type(fgumi_raw_bam::aux_data_slice(rec.as_ref()), tag).is_some()
        };
        // Forward, mapped read: 20M at 0-based pos 99 (alignment start 100).
        let mapped_fwd = || {
            fgumi_raw_bam::SamBuilder::new()
                .read_name(b"t")
                .sequence(b"ACGTACGTACGTACGTACGT")
                .qualities(&[30; 20])
                .cigar_ops(&[20u32 << 4])
                .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .build()
        };
        // Reverse, mapped read: 20M at 0-based pos 199 (alignment start 200).
        let mapped_rev = || {
            fgumi_raw_bam::SamBuilder::new()
                .read_name(b"t")
                .sequence(b"ACGTACGTACGTACGTACGT")
                .qualities(&[30; 20])
                .cigar_ops(&[20u32 << 4])
                .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
                .ref_id(0)
                .pos(199)
                .mapq(40)
                .build()
        };
        // Unmapped read as produced by makeReadUnmapped: no ref/pos, no strand, no CIGAR.
        let unmapped = |last: bool| {
            let seg = if last { raw_flags::LAST_SEGMENT } else { raw_flags::FIRST_SEGMENT };
            fgumi_raw_bam::SamBuilder::new()
                .read_name(b"t")
                .sequence(b"ACGTACGTACGTACGTACGT")
                .qualities(&[30; 20])
                .cigar_ops(&[])
                .flags(raw_flags::PAIRED | seg | raw_flags::UNMAPPED)
                .ref_id(-1)
                .pos(-1)
                .mapq(0)
                .build()
        };

        // --- Branch 1: both mapped ---
        let (mut r1, mut r2) = (mapped_fwd(), mapped_rev());
        set_mate_info_raw(&mut r1, &mut r2);
        assert_eq!((r1.mate_ref_id(), r1.mate_pos()), (0, 199), "r1 mate coords <- r2");
        assert!(has_flag(&r1, raw_flags::MATE_REVERSE), "r1 mate-reverse (r2 is reverse)");
        assert!(!has_flag(&r1, raw_flags::MATE_UNMAPPED), "r1 mate mapped");
        assert_eq!((r2.mate_ref_id(), r2.mate_pos()), (0, 99), "r2 mate coords <- r1");
        assert!(!has_flag(&r2, raw_flags::MATE_REVERSE), "r2 mate not reverse (r1 forward)");
        // TLEN: fwd 5' = 100, rev 5' = alignmentEnd 219 -> 219 - 100 + 1 = 120.
        assert_eq!(r1.template_length(), 120, "r1 TLEN");
        assert_eq!(r2.template_length(), -120, "r2 TLEN");
        assert!(has_tag(&r1, *SamTag::MQ) && has_tag(&r1, *SamTag::MC), "r1 MQ/MC set");
        assert!(has_tag(&r2, *SamTag::MQ) && has_tag(&r2, *SamTag::MC), "r2 MQ/MC set");

        // --- Branch 2: one mapped (r1), one unmapped (r2) ---
        let (mut r1, mut r2) = (mapped_fwd(), unmapped(true));
        set_mate_info_raw(&mut r1, &mut r2);
        // Unmapped r2 is relocated to r1's coordinate.
        assert_eq!((r2.ref_id(), r2.pos()), (0, 99), "unmapped r2 placed at mate coord");
        assert_eq!((r2.mate_ref_id(), r2.mate_pos()), (0, 99), "r2 mate coords <- r1");
        assert!(!has_flag(&r2, raw_flags::MATE_UNMAPPED), "r2's mate (r1) is mapped");
        assert!(
            has_tag(&r2, *SamTag::MQ) && has_tag(&r2, *SamTag::MC),
            "r2 MQ/MC from mapped mate"
        );
        // Mapped r1 points at the co-located unmapped mate.
        assert_eq!((r1.mate_ref_id(), r1.mate_pos()), (0, 99), "r1 mate coords <- unmapped r2");
        assert!(has_flag(&r1, raw_flags::MATE_UNMAPPED), "r1 mate unmapped");
        assert!(!has_flag(&r1, raw_flags::MATE_REVERSE), "r1 mate-reverse cleared (unmapped)");
        assert!(!has_tag(&r1, *SamTag::MQ) && !has_tag(&r1, *SamTag::MC), "r1 MQ/MC removed");
        assert_eq!(
            (r1.template_length(), r2.template_length()),
            (0, 0),
            "TLEN 0 when a mate unmapped"
        );

        // --- Branch 3: both unmapped ---
        let (mut r1, mut r2) = (unmapped(false), unmapped(true));
        set_mate_info_raw(&mut r1, &mut r2);
        for rec in [&r1, &r2] {
            assert_eq!((rec.ref_id(), rec.pos()), (-1, -1), "unmapped coords cleared");
            assert_eq!(
                (rec.mate_ref_id(), rec.mate_pos()),
                (-1, -1),
                "unmapped mate coords cleared"
            );
            assert!(has_flag(rec, raw_flags::MATE_UNMAPPED), "mate-unmapped set");
            assert!(!has_flag(rec, raw_flags::MATE_REVERSE), "mate-reverse cleared");
            assert!(!has_tag(rec, *SamTag::MQ) && !has_tag(rec, *SamTag::MC), "MQ/MC removed");
            assert_eq!(rec.template_length(), 0, "TLEN 0");
        }
    }

    /// A read that clipping left unmapped is relocated to its mapped mate's
    /// coordinate by `set_mate_info_raw` (htsjdk `SamPairUtil.setMateInfo`). fgbio
    /// emits the placed position's bin for such a read (htsjdk recomputes the
    /// indexing bin on write); the raw pipeline must do so explicitly.
    #[test]
    fn set_mate_info_raw_recomputes_bin_for_relocated_unmapped_read() {
        use fgumi_raw_bam::flags as raw_flags;
        let mut mapped = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"t")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(40)
            .build();
        let mut unmapped = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"t")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[])
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::UNMAPPED)
            .ref_id(-1)
            .pos(-1)
            .mapq(0)
            .build();
        assert_eq!(unmapped.bin(), fgumi_raw_bam::UNMAPPED_BIN, "precondition: unmapped bin");

        set_mate_info_raw(&mut mapped, &mut unmapped);

        assert_eq!(unmapped.pos(), 99, "precondition: relocated to mate coord");
        assert_eq!(
            unmapped.bin(),
            fgumi_raw_bam::reg2bin(99, 100),
            "relocated unmapped read must carry its placed position's bin, not the stale 4680"
        );
    }

    /// htsjdk's one-mapped branch sets the mapped read's mate-reverse flag from the
    /// unmapped read's *actual current* strand (`SamPairUtil.java:267`), not an assumed
    /// `false`. This only diverges for a read that was already unmapped on input while
    /// still carrying the REVERSE flag — the clipper's own unmap clears it, but
    /// `set_mate_info_raw` runs on every pair, so such a read reaches this branch.
    #[test]
    fn test_set_mate_info_raw_one_mapped_uses_unmapped_reads_strand() {
        use fgumi_raw_bam::flags as raw_flags;
        let has_flag = |rec: &RawRecord, f: u16| rec.flags() & f != 0;

        // Mapped forward r1.
        let mut r1 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"t")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[20u32 << 4])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .build();
        // Unmapped r2 that still carries the REVERSE flag (already unmapped on input).
        let mut r2 = fgumi_raw_bam::SamBuilder::new()
            .read_name(b"t")
            .sequence(b"ACGTACGTACGTACGTACGT")
            .qualities(&[30; 20])
            .cigar_ops(&[])
            .flags(
                raw_flags::PAIRED
                    | raw_flags::LAST_SEGMENT
                    | raw_flags::UNMAPPED
                    | raw_flags::REVERSE,
            )
            .ref_id(-1)
            .pos(-1)
            .mapq(0)
            .build();

        set_mate_info_raw(&mut r1, &mut r2);

        assert!(
            has_flag(&r1, raw_flags::MATE_REVERSE),
            "mapped r1 mate-reverse must reflect unmapped r2's actual REVERSE flag"
        );
    }

    /// Ports of fgbio `ClipBamTest` (fgbio commit `e51a661`) that fgumi's own `clip` tests above
    /// either did not cover or covered more weakly (different inputs, missing assertions).
    /// Inputs and expected values are copied verbatim from the fgbio tests; each test cites its
    /// source as `ClipBamTest.scala:<line>`, and any fixture deviation is explained there.
    ///
    /// fgbio's `clipper.clipPair(r1, r2)` cases run [`ClipParams::clip_pair`] with a `Hard`
    /// clipper (fgbio `ClipBam`'s default mode); its `.execute()` cases run [`Clip::execute`]
    /// end to end.
    mod fgbio_clip_bam_tests {
        use super::*;
        use crate::metrics::clip::{ClippingMetrics, ReadType};
        use crate::sam::builder::PairBuilder;
        use fgumi_raw_bam::{encode_record_buf_to_raw, raw_record_to_record_buf};
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record::cigar::op::Kind;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value;

        /// Length of fgbio's `chr1` test reference, which is all `A`.
        const REFERENCE_LENGTH: usize = 5000;

        /// Writes fgbio's test reference (`chr1`: 5000 `A`s) plus its `.fai`; `NM`/`UQ`/`MD`
        /// expectations in these tests are computed against it.
        fn write_all_a_reference(dir: &TempDir) -> PathBuf {
            let path = dir.path().join("ref.fa");
            std::fs::write(&path, format!(">chr1\n{}\n", "A".repeat(REFERENCE_LENGTH)))
                .expect("write reference");
            std::fs::write(
                dir.path().join("ref.fa.fai"),
                format!(
                    "chr1\t{REFERENCE_LENGTH}\t6\t{REFERENCE_LENGTH}\t{}\n",
                    REFERENCE_LENGTH + 1
                ),
            )
            .expect("write reference index");
            path
        }

        /// fgbio `new SamBuilder(readLength).addPair(...)`: builds an FR pair (by default) of
        /// all-`A` reads of `read_length` bases, letting `configure` set positions, strands and
        /// CIGARs, and encodes both reads as [`RawRecord`]s.
        fn fgbio_pair(
            read_length: usize,
            configure: impl for<'a> FnOnce(PairBuilder<'a>) -> PairBuilder<'a>,
        ) -> (RawRecord, RawRecord) {
            let mut builder = SamBuilder::new();
            let bases = "A".repeat(read_length);
            let pair = builder.add_pair().name("q").bases1(&bases).bases2(&bases);
            let (r1, r2) = configure(pair).build();
            let header = builder.header.clone();
            (
                encode_record_buf_to_raw(&r1, &header).expect("encode r1"),
                encode_record_buf_to_raw(&r2, &header).expect("encode r2"),
            )
        }

        /// fgbio `clipper.clipPair(r1, r2)`: runs the per-template clipping that `clip`'s
        /// options configure on one pair, in fgbio `ClipBam`'s default `Hard` mode.
        fn clip_pair_hard(clip: &Clip, r1: &mut RawRecord, r2: &mut RawRecord) {
            ClipParams::from_clip(clip)
                .clip_pair(&RawRecordClipper::new(ClippingMode::Hard), r1, r2, None)
                .expect("clip_pair should succeed");
        }

        /// A `Clip` with only `--clip-overlapping-reads` (and no fixed clipping) enabled.
        fn overlap_only() -> Clip {
            let mut clip = make_clip(0, 0, 0, 0);
            clip.clip_overlapping_reads = true;
            clip
        }

        /// 1-based alignment start of a mapped raw record, as fgbio's `rec.start`.
        fn start(rec: &RawRecord) -> usize {
            rec.alignment_start_1based().expect("record should be mapped")
        }

        /// 1-based inclusive alignment end of a mapped raw record, as fgbio's `rec.end`.
        fn end(rec: &RawRecord) -> usize {
            rec.alignment_end_1based().expect("record should be mapped")
        }

        /// fgbio `StartAndEnd.checkClipping`: asserts `rec` lost exactly `five_prime` aligned
        /// bases at its 5' end and `three_prime` at its 3' end relative to `prior` `(start,
        /// end)`, accounting for strand.
        fn assert_clipped_by(
            prior: (usize, usize),
            rec: &RawRecord,
            five_prime: usize,
            three_prime: usize,
            label: &str,
        ) {
            let (prior_start, prior_end) = prior;
            let expected = if rec.is_reverse() {
                (prior_start + three_prime, prior_end - five_prime)
            } else {
                (prior_start + five_prime, prior_end - three_prime)
            };
            assert_eq!((start(rec), end(rec)), expected, "{label} (start, end)");
        }

        /// `ClipBamTest.scala:86` "not clip reads where either read is unaligned".
        #[test]
        fn clip_pair_does_not_clip_when_a_read_is_unaligned() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| p.start1(100).unmapped2());
            let expected = r1.cigar_to_string();
            clip_pair_hard(&overlap_only(), &mut r1, &mut r2);
            assert_eq!(r1.cigar_to_string(), expected);
        }

        /// `ClipBamTest.scala:95` "not clip reads that are on different chromosomes".
        #[test]
        fn clip_pair_does_not_clip_reads_on_different_chromosomes() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| p.start1(100).start2(100).contig2(1));
            let expected = (r1.cigar_to_string(), r2.cigar_to_string());
            clip_pair_hard(&overlap_only(), &mut r1, &mut r2);
            assert_eq!((r1.cigar_to_string(), r2.cigar_to_string()), expected);
        }

        /// `ClipBamTest.scala:108` "not clip reads that are abutting but not overlapped".
        #[test]
        fn clip_pair_does_not_clip_abutting_reads() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| p.start1(100).start2(150));
            let expected = (r1.cigar_to_string(), r2.cigar_to_string());
            clip_pair_hard(&overlap_only(), &mut r1, &mut r2);
            assert_eq!((r1.cigar_to_string(), r2.cigar_to_string()), expected);
        }

        /// `ClipBamTest.scala:119` "not clip non-FR reads".
        #[test]
        fn clip_pair_does_not_clip_non_fr_reads() {
            let (mut r1, mut r2) =
                fgbio_pair(50, |p| p.start1(100).start2(100).strand2(Strand::Plus));
            let expected = (r1.cigar_to_string(), r2.cigar_to_string());
            clip_pair_hard(&overlap_only(), &mut r1, &mut r2);
            assert_eq!((r1.cigar_to_string(), r2.cigar_to_string()), expected);
        }

        /// `ClipBamTest.scala:130` "clip reads that are fully overlapped".
        #[test]
        fn clip_pair_hard_clips_fully_overlapped_reads() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| p.start1(100).start2(100));
            clip_pair_hard(&overlap_only(), &mut r1, &mut r2);
            assert_eq!(r1.cigar_to_string(), "25M25H");
            assert_eq!(r2.cigar_to_string(), "25H25M");
        }

        /// `ClipBamTest.scala:160` "handle reads that contain deletions".
        #[test]
        fn clip_pair_removes_overlap_of_reads_with_deletions() {
            let (mut r1, mut r2) =
                fgbio_pair(50, |p| p.start1(100).start2(130).cigar1("40M2D10M").cigar2("10M2D40M"));
            assert!(end(&r1) >= start(&r2), "fixture should overlap");
            clip_pair_hard(&overlap_only(), &mut r1, &mut r2);
            assert!(end(&r1) < start(&r2), "r1.end={} r2.start={}", end(&r1), start(&r2));
        }

        /// `ClipBamTest.scala:170` "clip a fixed amount on the ends of the reads with reads that
        /// do not overlap".
        #[test]
        fn clip_pair_fixed_clipping_without_overlap() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| p.start1(100).start2(150));
            let (prior1, prior2) = ((start(&r1), end(&r1)), (start(&r2), end(&r2)));
            assert_eq!(end(&r1), start(&r2) - 1);
            clip_pair_hard(&make_clip(1, 2, 3, 4), &mut r1, &mut r2);
            assert_clipped_by(prior1, &r1, 1, 2, "r1");
            assert_clipped_by(prior2, &r2, 3, 4, "r2");
        }

        /// `ClipBamTest.scala:182` "clip a fixed amount on the ends of the reads with reads with
        /// clipping present": existing hard clipping counts toward the fixed amounts.
        ///
        /// Deviation: fgbio sets `4H46M` / `44M6H` on 50-base reads, so its SEQ is longer than
        /// the CIGAR's query length. fgumi encodes records to raw BAM, which rejects that, so
        /// the reads here carry 46 / 44 bases to match their CIGARs. The expected clipping is
        /// unchanged.
        #[test]
        fn clip_pair_fixed_clipping_counts_existing_hard_clips() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| {
                p.start1(104)
                    .start2(150)
                    .cigar1("4H46M")
                    .bases1(&"A".repeat(46))
                    .cigar2("44M6H")
                    .bases2(&"A".repeat(44))
            });
            let (prior1, prior2) = ((start(&r1), end(&r1)), (start(&r2), end(&r2)));
            assert_eq!(end(&r1), start(&r2) - 1);
            clip_pair_hard(&make_clip(5, 2, 3, 4), &mut r1, &mut r2);
            // R1: one more 5' base (4H already counts toward 5), two 3' bases.
            assert_clipped_by(prior1, &r1, 1, 2, "r1");
            // R2: no more 5' bases (6H already exceeds 3), four 3' bases.
            assert_clipped_by(prior2, &r2, 0, 4, "r2");
        }

        /// `ClipBamTest.scala:201` "clip a fixed amount on the ends of the reads then clip
        /// overlapping reads": fixed clipping is applied first, then the remaining overlap.
        #[test]
        fn clip_pair_fixed_clipping_then_overlap_clipping() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| p.start1(100).start2(146));
            let (prior1, prior2) = ((start(&r1), end(&r1)), (start(&r2), end(&r2)));
            assert_eq!(end(&r1), start(&r2) + 3, "four bases overlap");
            let mut clip = make_clip(0, 1, 0, 1);
            clip.clip_overlapping_reads = true;
            clip_pair_hard(&clip, &mut r1, &mut r2);
            assert_eq!(end(&r1), start(&r2) - 1);
            assert_clipped_by(prior1, &r1, 0, 2, "r1");
            assert_clipped_by(prior2, &r2, 0, 2, "r2");
        }

        /// `ClipBamTest.scala:216` "clip a fixed amount on the ends of the reads in
        /// $strand1/$strand2", for all four strand combinations.
        #[rstest]
        #[case::plus_plus(Strand::Plus, Strand::Plus)]
        #[case::plus_minus(Strand::Plus, Strand::Minus)]
        #[case::minus_plus(Strand::Minus, Strand::Plus)]
        #[case::minus_minus(Strand::Minus, Strand::Minus)]
        fn clip_pair_fixed_clipping_is_strand_aware(
            #[case] strand1: Strand,
            #[case] strand2: Strand,
        ) {
            let (mut r1, mut r2) =
                fgbio_pair(50, |p| p.start1(100).start2(150).strand1(strand1).strand2(strand2));
            let (prior1, prior2) = ((start(&r1), end(&r1)), (start(&r2), end(&r2)));
            assert_eq!(end(&r1), start(&r2) - 1);
            clip_pair_hard(&make_clip(1, 2, 3, 4), &mut r1, &mut r2);
            assert_clipped_by(prior1, &r1, 1, 2, "r1");
            assert_clipped_by(prior2, &r2, 3, 4, "r2");
        }

        /// All 15 counters of a [`ClippingMetrics`], in fgbio's column order.
        fn metric_values(m: &ClippingMetrics) -> [usize; 15] {
            [
                m.reads,
                m.reads_unmapped,
                m.reads_clipped_pre,
                m.reads_clipped_post,
                m.reads_clipped_five_prime,
                m.reads_clipped_three_prime,
                m.reads_clipped_overlapping,
                m.reads_clipped_extending,
                m.bases,
                m.bases_clipped_pre,
                m.bases_clipped_post,
                m.bases_clipped_five_prime,
                m.bases_clipped_three_prime,
                m.bases_clipped_overlapping,
                m.bases_clipped_extending,
            ]
        }

        /// A [`ClippingMetrics`] whose 15 counters are `first, first + 1, ..., first + 14`, as
        /// in fgbio's `ClippingMetrics(readType, 1, 2, 3, ...)` fixtures.
        fn metrics_counting_up_from(read_type: ReadType, first: usize) -> ClippingMetrics {
            let mut m = ClippingMetrics::new(read_type);
            let fields = [
                &mut m.reads,
                &mut m.reads_unmapped,
                &mut m.reads_clipped_pre,
                &mut m.reads_clipped_post,
                &mut m.reads_clipped_five_prime,
                &mut m.reads_clipped_three_prime,
                &mut m.reads_clipped_overlapping,
                &mut m.reads_clipped_extending,
                &mut m.bases,
                &mut m.bases_clipped_pre,
                &mut m.bases_clipped_post,
                &mut m.bases_clipped_five_prime,
                &mut m.bases_clipped_three_prime,
                &mut m.bases_clipped_overlapping,
                &mut m.bases_clipped_extending,
            ];
            for (offset, field) in fields.into_iter().enumerate() {
                *field = first + offset;
            }
            m
        }

        /// `ClipBamTest.scala:230` "add two metrics" (`ClippingMetrics.add`): every counter is
        /// summed, and the sum differs from each input in every field. fgbio builds `readOne`
        /// with `ReadType.ReadTwo`; that is kept.
        #[test]
        fn clipping_metrics_add_sums_every_field() {
            let fragment = metrics_counting_up_from(ReadType::Fragment, 1);
            let read_one = metrics_counting_up_from(ReadType::ReadTwo, 2);
            let read_two = metrics_counting_up_from(ReadType::ReadTwo, 3);
            let mut added = ClippingMetrics::new(ReadType::All);
            added.add(&fragment);
            added.add(&read_one);
            added.add(&read_two);

            assert_eq!(added.read_type, ReadType::All);
            assert_eq!(
                metric_values(&added),
                [6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48]
            );
            assert_ne!(added.read_type, fragment.read_type);
            for (field, (sum, input)) in
                metric_values(&added).into_iter().zip(metric_values(&fragment)).enumerate()
            {
                assert_ne!(sum, input, "field {field} should differ from the fragment input");
            }
        }

        /// The `(NM, UQ, MD)` htsjdk's `calculateMdAndNmTags` / `sumQualitiesOfMismatches` give
        /// `rec` against fgbio's all-`A` reference. Supports only the CIGAR operators these
        /// fixtures use.
        fn expected_nm_uq_md(rec: &RecordBuf) -> (i64, i64, String) {
            use std::fmt::Write as _;
            let bases: &[u8] = rec.sequence().as_ref();
            let quals: &[u8] = rec.quality_scores().as_ref();
            let (mut offset, mut nm, mut uq, mut matches) = (0usize, 0i64, 0i64, 0usize);
            let mut md = String::new();
            for op in rec.cigar().as_ref() {
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        for _ in 0..op.len() {
                            if bases[offset].eq_ignore_ascii_case(&b'A') {
                                matches += 1;
                            } else {
                                nm += 1;
                                uq += i64::from(quals[offset]);
                                write!(md, "{matches}A").expect("write to String");
                                matches = 0;
                            }
                            offset += 1;
                        }
                    }
                    Kind::Insertion | Kind::SoftClip => offset += op.len(),
                    Kind::HardClip => {}
                    kind => panic!("unsupported CIGAR operator in fixture: {kind:?}"),
                }
            }
            write!(md, "{matches}").expect("write to String");
            (nm, uq, md)
        }

        /// Sets `NM`, `UQ` and `MD` on `rec` to their correct values, as fgbio's fixtures do
        /// before running `ClipBam`.
        fn set_nm_uq_md(rec: &mut RecordBuf) {
            let (nm, uq, md) = expected_nm_uq_md(rec);
            let data = rec.data_mut();
            data.insert(
                Tag::from(SamTag::NM),
                Value::from(i32::try_from(nm).expect("NM fits i32")),
            );
            data.insert(
                Tag::from(SamTag::UQ),
                Value::from(i32::try_from(uq).expect("UQ fits i32")),
            );
            data.insert(Tag::from(SamTag::MD), Value::from(md.as_str()));
        }

        /// Asserts `rec` carries `NM`, `UQ` and `MD` equal to [`expected_nm_uq_md`].
        fn assert_nm_uq_md_recomputed(rec: &RecordBuf) {
            let (nm, uq, md) = expected_nm_uq_md(rec);
            let int_tag = |tag: SamTag| {
                rec.data()
                    .get(&Tag::from(tag))
                    .and_then(Value::as_int)
                    .unwrap_or_else(|| panic!("{tag:?} missing or not an integer"))
            };
            let md_tag = match rec.data().get(&Tag::from(SamTag::MD)) {
                Some(Value::String(s)) => s.to_string(),
                other => panic!("MD missing or not a string: {other:?}"),
            };
            let label = cigar_string(rec);
            assert_eq!((int_tag(SamTag::NM), int_tag(SamTag::UQ), md_tag), (nm, uq, md), "{label}");
        }

        /// The 50-base all-`A` read fgbio's NM/UQ/MD fixtures use, with a single `C` at
        /// `mismatch_at`.
        fn all_a_with_mismatch(mismatch_at: usize) -> String {
            let mut bases = vec![b'A'; 50];
            bases[mismatch_at] = b'C';
            String::from_utf8(bases).expect("ASCII bases")
        }

        /// The first 12 values of Java's `new Random(1).nextInt(50)`: the mismatch positions
        /// fgbio's NM/UQ/MD fixtures draw, in the order its `SamBuilder` yields the reads.
        const FGBIO_MISMATCH_POSITIONS: [usize; 12] = [35, 38, 47, 13, 4, 4, 34, 6, 28, 48, 19, 23];

        /// A `Clip` reading `input` and writing `output` (with `metrics`) against `reference`,
        /// in `mode`, with all clipping disabled; callers enable what their case needs.
        fn end_to_end_clip(
            input: PathBuf,
            output: PathBuf,
            reference: PathBuf,
            metrics: PathBuf,
            mode: ClippingMode,
        ) -> Clip {
            let mut clip = make_clip(0, 0, 0, 0);
            clip.io.input = input;
            clip.io.output = output;
            clip.reference = reference;
            clip.metrics = Some(metrics);
            clip.clipping_mode = mode;
            clip
        }

        /// Reads `clip`'s metrics file back.
        fn read_clipping_metrics(path: &std::path::Path) -> Vec<ClippingMetrics> {
            crate::metrics::read_metrics(path, "clipping").expect("read clipping metrics")
        }

        /// The first and last CIGAR operations of a mapped record.
        fn first_and_last_ops(rec: &RecordBuf) -> (Op, Op) {
            let ops = rec.cigar().as_ref();
            (*ops.first().expect("non-empty CIGAR"), *ops.last().expect("non-empty CIGAR"))
        }

        /// `ClipBamTest.scala:245` "clip overlapping reads, update mate info, and reset NM, UQ &
        /// MD": six pairs overlapping by 10, 8, 6, 4, 2 and 0 bases, with fixed 5' clipping (R1
        /// 2, R2 3) and overlap clipping in `Hard` mode. Checks clipping, coordinate order, mate
        /// info (`MC`), recomputed `NM`/`UQ`/`MD`, and every metrics row.
        #[test]
        fn clip_execute_pairs_updates_mate_info_tags_and_metrics() {
            let dir = TempDir::new().expect("temp dir");
            let reference = write_all_a_reference(&dir);
            let mut source = SamBuilder::with_single_ref("chr1", REFERENCE_LENGTH);
            let mut input = SamBuilder::with_single_ref("chr1", REFERENCE_LENGTH);
            input.set_queryname_sort_order();
            let starts = [(100, 140), (200, 242), (300, 344), (400, 446), (500, 548), (600, 650)];
            for (i, (start1, start2)) in starts.into_iter().enumerate() {
                let (mut r1, mut r2) = source
                    .add_pair()
                    .name(&format!("q{}", i + 1))
                    .bases1(&all_a_with_mismatch(FGBIO_MISMATCH_POSITIONS[2 * i]))
                    .bases2(&all_a_with_mismatch(FGBIO_MISMATCH_POSITIONS[2 * i + 1]))
                    .start1(start1)
                    .start2(start2)
                    .build();
                set_nm_uq_md(&mut r1);
                set_nm_uq_md(&mut r2);
                input.push_record(r1);
                input.push_record(r2);
            }
            let (input_path, output, metrics) =
                (dir.path().join("in.bam"), dir.path().join("out.bam"), dir.path().join("m.txt"));
            input.write(&input_path).expect("write input");
            let mut clip = end_to_end_clip(
                input_path,
                output.clone(),
                reference,
                metrics.clone(),
                ClippingMode::Hard,
            );
            clip.read_one_five_prime = 2;
            clip.read_two_five_prime = 3;
            clip.clip_overlapping_reads = true;
            clip.execute("test").expect("clip should succeed");

            let clipped = read_bam_records(&output).expect("read output");
            assert_eq!(clipped.len(), 12);
            let is_hard = |op: Op| op.kind() == Kind::HardClip;
            let is_hard_of = |op: Op, len: usize| is_hard(op) && op.len() == len;
            let (minus, plus): (Vec<_>, Vec<_>) =
                clipped.iter().partition(|r| r.flags().is_reverse_complemented());
            let count = |recs: &[&RecordBuf], pred: &dyn Fn((Op, Op)) -> bool| {
                recs.iter().filter(|r| pred(first_and_last_ops(r))).count()
            };
            // Overlap clipping hit every pair but q6.
            assert_eq!(count(&minus, &|(first, _)| is_hard(first)), 5);
            assert_eq!(count(&plus, &|(_, last)| is_hard(last)), 5);
            // Fixed 5' clipping hit every read.
            assert_eq!(count(&minus, &|(_, last)| is_hard_of(last, 3)), 6);
            assert_eq!(count(&plus, &|(first, _)| is_hard_of(first, 2)), 6);

            for pair in clipped.windows(2) {
                assert!(pair[0].alignment_start() <= pair[1].alignment_start(), "coordinate order");
            }
            let mut templates: std::collections::BTreeMap<_, Vec<&RecordBuf>> =
                std::collections::BTreeMap::new();
            for rec in &clipped {
                templates.entry(rec.name()).or_default().push(rec);
            }
            assert_eq!(templates.len(), 6);
            for template in templates.values() {
                let [lhs, rhs] = template.as_slice() else {
                    panic!("expected a pair, got {template:?}")
                };
                for (a, b) in [(lhs, rhs), (rhs, lhs)] {
                    assert_eq!(a.mate_alignment_start(), b.alignment_start(), "mate start");
                    assert_eq!(
                        a.flags().is_mate_reverse_complemented(),
                        b.flags().is_reverse_complemented(),
                        "mate strand"
                    );
                    match a.data().get(&Tag::from(SamTag::MC)) {
                        Some(Value::String(mc)) => assert_eq!(mc.to_string(), cigar_string(b)),
                        other => panic!("MC missing or not a string: {other:?}"),
                    }
                }
            }
            for rec in &clipped {
                assert_nm_uq_md_recomputed(rec);
            }

            let rows = read_clipping_metrics(&metrics);
            let row = |read_type: ReadType| {
                let row = rows.iter().find(|m| m.read_type == read_type);
                metric_values(row.unwrap_or_else(|| panic!("no {read_type:?} row")))
            };
            // Columns: reads, unmapped, clipped pre/post/5'/3'/overlapping/extending, bases,
            // bases clipped pre/post/5'/3'/overlapping/extending.
            let read_one = [6, 0, 0, 6, 6, 0, 5, 0, 273, 0, 27, 12, 0, 15, 0];
            let read_two = [6, 0, 0, 6, 6, 0, 5, 0, 267, 0, 33, 18, 0, 15, 0];
            let pair: [usize; 15] = std::array::from_fn(|i| read_one[i] + read_two[i]);
            assert_eq!(rows.len(), 5);
            assert_eq!(row(ReadType::Fragment), [0; 15], "Fragment");
            assert_eq!(row(ReadType::ReadOne), read_one, "ReadOne");
            assert_eq!(row(ReadType::ReadTwo), read_two, "ReadTwo");
            assert_eq!(row(ReadType::Pair), pair, "Pair");
            assert_eq!(row(ReadType::All), pair, "All");
        }

        /// `ClipBamTest.scala:351` "clip fragment reads, and reset NM, UQ & MD": three fragments
        /// with R1 fixed clipping (5' 2, 3' 10); R2 clipping and overlap clipping do not apply to
        /// fragments, and the `40S10M` fragment is clipped away entirely and unmapped.
        #[test]
        fn clip_execute_fragments_resets_tags_and_metrics() {
            let dir = TempDir::new().expect("temp dir");
            let reference = write_all_a_reference(&dir);
            let mut source = SamBuilder::with_single_ref("chr1", REFERENCE_LENGTH);
            let mut input = SamBuilder::with_single_ref("chr1", REFERENCE_LENGTH);
            input.set_queryname_sort_order();
            let fragments = [
                (100, Strand::Plus, "50M"),
                (200, Strand::Minus, "50M"),
                (300, Strand::Plus, "40S10M"),
            ];
            for (i, (start, strand, cigar)) in fragments.into_iter().enumerate() {
                let mut rec = source
                    .add_frag()
                    .name(&format!("f{}", i + 1))
                    .bases(&all_a_with_mismatch(FGBIO_MISMATCH_POSITIONS[i]))
                    .start(start)
                    .strand(strand)
                    .cigar(cigar)
                    .build();
                set_nm_uq_md(&mut rec);
                input.push_record(rec);
            }
            let (input_path, output, metrics) =
                (dir.path().join("in.bam"), dir.path().join("out.bam"), dir.path().join("m.txt"));
            input.write(&input_path).expect("write input");
            let mut clip = end_to_end_clip(
                input_path,
                output.clone(),
                reference,
                metrics.clone(),
                ClippingMode::Hard,
            );
            (clip.read_one_five_prime, clip.read_one_three_prime) = (2, 10);
            (clip.read_two_five_prime, clip.read_two_three_prime) = (5, 5);
            clip.clip_overlapping_reads = true;
            clip.execute("test").expect("clip should succeed");

            let clipped = read_bam_records(&output).expect("read output");
            assert_eq!(clipped.len(), 3);
            let mapped: Vec<_> = clipped.iter().filter(|r| !r.flags().is_unmapped()).collect();
            let is_2h = |op: Op| op.kind() == Kind::HardClip && op.len() == 2;
            let minus_5p_clipped = mapped
                .iter()
                .filter(|r| r.flags().is_reverse_complemented() && is_2h(first_and_last_ops(r).1))
                .count();
            let plus_5p_clipped = mapped
                .iter()
                .filter(|r| !r.flags().is_reverse_complemented() && is_2h(first_and_last_ops(r).0))
                .count();
            assert_eq!((minus_5p_clipped, plus_5p_clipped), (1, 1));
            for pair in clipped.windows(2) {
                let (lhs, rhs) = (&pair[0], &pair[1]);
                if !lhs.flags().is_unmapped() && !rhs.flags().is_unmapped() {
                    assert!(lhs.alignment_start() <= rhs.alignment_start(), "coordinate order");
                } else if lhs.flags().is_unmapped() {
                    assert!(rhs.flags().is_unmapped(), "unmapped reads sort last");
                }
            }
            for rec in &mapped {
                assert_nm_uq_md_recomputed(rec);
            }

            let rows = read_clipping_metrics(&metrics);
            // Columns as in `clip_execute_pairs_updates_mate_info_tags_and_metrics`.
            let fragment = [3, 1, 1, 3, 2, 3, 0, 0, 76, 40, 74, 4, 30, 0, 0];
            assert_eq!(rows.len(), 5);
            for row in &rows {
                let expected = match row.read_type {
                    ReadType::Fragment | ReadType::All => fragment,
                    ReadType::ReadOne | ReadType::ReadTwo | ReadType::Pair => [0; 15],
                };
                assert_eq!(metric_values(row), expected, "{:?}", row.read_type);
            }
        }

        /// fgbio's `--upgrade-clipping` fixture: two 50 bp fragments, `q1` (`+`, start 100) and
        /// `q2` (`-`, start 200), each clipped 10 bases at the 5' end and 4 at the 3' end in
        /// `prior` mode (with auto-clip attributes). With `reclip_in_prior_mode`, each is then
        /// clipped again (5' 5, 3' 2), which must be a no-op, as in fgbio's upgrade cases.
        /// Returns the fragments as written to `clip`'s input, and `clip`'s output in `mode`.
        fn run_upgrade_clipping(
            prior: ClippingMode,
            mode: ClippingMode,
            reclip_in_prior_mode: bool,
        ) -> (Vec<RecordBuf>, Vec<RecordBuf>) {
            let dir = TempDir::new().expect("temp dir");
            let reference = write_all_a_reference(&dir);
            let mut source = SamBuilder::with_single_ref("chr1", REFERENCE_LENGTH);
            let mut input = SamBuilder::with_single_ref("chr1", REFERENCE_LENGTH);
            input.set_queryname_sort_order();
            let quals: Vec<u8> = (0..50u8).map(|i| 10 + i % 30).collect();
            let fragments = [
                ("q1", "ACGTTGCAAC".repeat(5), 100, Strand::Plus),
                ("q2", "TTGACCAGTA".repeat(5), 200, Strand::Minus),
            ];
            let header = source.header.clone();
            let clipper = RawRecordClipper::with_auto_clip(prior, true);
            for (name, bases, start, strand) in fragments {
                let frag = source
                    .add_frag()
                    .name(name)
                    .bases(&bases)
                    .quals(&quals)
                    .start(start)
                    .strand(strand)
                    .build();
                let mut raw = encode_record_buf_to_raw(&frag, &header).expect("encode");
                assert_eq!(clipper.clip_5_prime_end_of_read_raw(&mut raw, 10), 10);
                assert_eq!(clipper.clip_3_prime_end_of_read_raw(&mut raw, 4), 4);
                if reclip_in_prior_mode {
                    assert_eq!(clipper.clip_5_prime_end_of_read_raw(&mut raw, 5), 0);
                    assert_eq!(clipper.clip_3_prime_end_of_read_raw(&mut raw, 2), 0);
                }
                input.push_record(raw_record_to_record_buf(&raw, &header).expect("decode"));
            }
            let (input_path, output) = (dir.path().join("in.bam"), dir.path().join("out.bam"));
            input.write(&input_path).expect("write input");
            let mut clip = end_to_end_clip(
                input_path,
                output.clone(),
                reference,
                dir.path().join("m.txt"),
                mode,
            );
            clip.upgrade_clipping = true;
            clip.execute("test").expect("clip should succeed");
            let clipped = read_bam_records(&output).expect("read output");
            assert_eq!(clipped.len(), 2);
            (input.records().to_vec(), clipped)
        }

        /// fgbio's `maskBases`/`maskQuals`: `values` with its first `leading` and last
        /// `trailing` entries replaced by `fill`.
        fn masked(values: &[u8], leading: usize, trailing: usize, fill: u8) -> Vec<u8> {
            let mut out = values.to_vec();
            let len = out.len();
            out[..leading].fill(fill);
            out[len - trailing..].fill(fill);
            out
        }

        /// Asserts `clipped` is `prior` with its clipping in `expected` mode, per fgbio's
        /// upgrade expectations: `q1` (`+`) carries 10 5' / 4 3' clipped bases (`10?36M4?`),
        /// `q2` (`-`) the mirror (`4?36M10?`). `prior_mode` says whether `prior`'s bases were
        /// already hard-clipped away.
        fn assert_clipping_in_mode(
            prior_mode: ClippingMode,
            expected: ClippingMode,
            prior: &[RecordBuf],
            clipped: &[RecordBuf],
        ) {
            let seq = |r: &RecordBuf| r.sequence().as_ref().to_vec();
            let qual = |r: &RecordBuf| r.quality_scores().as_ref().to_vec();
            for (prior, clipped, (leading, trailing)) in
                [(&prior[0], &clipped[0], (10, 4)), (&prior[1], &clipped[1], (4, 10))]
            {
                let op = if expected == ClippingMode::Hard { 'H' } else { 'S' };
                assert_eq!(cigar_string(clipped), format!("{leading}{op}36M{trailing}{op}"));
                let (want_seq, want_qual) = match expected {
                    ClippingMode::Soft => (seq(prior), qual(prior)),
                    ClippingMode::SoftWithMask => (
                        masked(&seq(prior), leading, trailing, b'N'),
                        masked(&qual(prior), leading, trailing, fgumi_dna::MIN_PHRED),
                    ),
                    ClippingMode::Hard if prior_mode == ClippingMode::Hard => {
                        (seq(prior), qual(prior))
                    }
                    ClippingMode::Hard => (
                        seq(prior)[leading..50 - trailing].to_vec(),
                        qual(prior)[leading..50 - trailing].to_vec(),
                    ),
                };
                assert_eq!(
                    seq(clipped),
                    want_seq,
                    "bases of {cigar}",
                    cigar = cigar_string(clipped)
                );
                assert_eq!(
                    qual(clipped),
                    want_qual,
                    "quals of {cigar}",
                    cigar = cigar_string(clipped)
                );
            }
        }

        /// `ClipBamTest.scala:428` "upgrade existing clipping from $prior to $mode with
        /// --upgrade-clipping".
        #[rstest]
        #[case::soft_to_soft_with_mask(ClippingMode::Soft, ClippingMode::SoftWithMask)]
        #[case::soft_to_hard(ClippingMode::Soft, ClippingMode::Hard)]
        #[case::soft_with_mask_to_hard(ClippingMode::SoftWithMask, ClippingMode::Hard)]
        fn clip_execute_upgrades_existing_clipping(
            #[case] prior: ClippingMode,
            #[case] mode: ClippingMode,
        ) {
            let (input, clipped) = run_upgrade_clipping(prior, mode, true);
            assert_clipping_in_mode(prior, mode, &input, &clipped);
        }

        /// `ClipBamTest.scala:473` "not upgrade existing clipping from $prior to $mode with
        /// --upgrade-clipping": clipping already at or above `mode` is left as it was.
        #[rstest]
        #[case::soft_to_soft(ClippingMode::Soft, ClippingMode::Soft)]
        #[case::soft_with_mask_to_soft(ClippingMode::SoftWithMask, ClippingMode::Soft)]
        #[case::soft_with_mask_to_soft_with_mask(
            ClippingMode::SoftWithMask,
            ClippingMode::SoftWithMask
        )]
        #[case::hard_to_hard(ClippingMode::Hard, ClippingMode::Hard)]
        #[case::hard_to_soft_with_mask(ClippingMode::Hard, ClippingMode::SoftWithMask)]
        #[case::hard_to_soft(ClippingMode::Hard, ClippingMode::Soft)]
        fn clip_execute_does_not_downgrade_existing_clipping(
            #[case] prior: ClippingMode,
            #[case] mode: ClippingMode,
        ) {
            let (input, clipped) = run_upgrade_clipping(prior, mode, false);
            assert_clipping_in_mode(prior, prior, &input, &clipped);
        }

        /// `ClipBamTest.scala:518` "clip FR reads that extend past the mate".
        #[test]
        fn clip_pair_clips_reads_extending_past_mate() {
            let (mut r1, mut r2) = fgbio_pair(50, |p| p.start1(100).start2(90));
            assert_eq!((end(&r1), end(&r2)), (149, 139));
            let mut clip = make_clip(0, 0, 0, 0);
            clip.clip_extending_past_mate = true;
            clip_pair_hard(&clip, &mut r1, &mut r2);
            assert_eq!(start(&r1), start(&r2));
            assert_eq!(end(&r1), end(&r2));
        }

        /// A `Clip` with fixed clipping `(r1 5', r1 3', r2 5', r2 3')`, past-mate clipping, and
        /// optionally overlap clipping.
        fn past_mate_clip(fixed: [usize; 4], clip_overlapping: bool) -> Clip {
            let mut clip = make_clip(fixed[0], fixed[1], fixed[2], fixed[3]);
            clip.clip_extending_past_mate = true;
            clip.clip_overlapping_reads = clip_overlapping;
            clip
        }

        /// fgbio's past-mate cases: an FR pair of `read_length` reads at `starts`, whose prior
        /// ends must be `prior_ends`, clipped by `clip`, giving `(r1 start, r1 end, r2 start, r2
        /// end)`.
        #[rstest]
        // ClipBamTest.scala:530 "clip FR reads that extend past their mate and remove overlap"
        #[case::past_mate_and_overlap(100, (100, 90), (199, 189), past_mate_clip([0, 0, 0, 0], true), (100, 144, 145, 189))]
        // ClipBamTest.scala:551 "clip FR reads that extend past their mate with asymmetrical five prime hard clipping"
        #[case::asymmetric_five_prime(200, (100, 90), (299, 289), past_mate_clip([10, 0, 50, 0], false), (110, 239, 110, 239))]
        // ClipBamTest.scala:574 "clip FR reads that extend past their mate with some irrelevant three prime clipping and removal of overlap"
        #[case::three_prime_and_overlap(200, (100, 90), (299, 289), past_mate_clip([0, 0, 0, 50], true), (100, 194, 195, 289))]
        // ClipBamTest.scala:597 "clip FR reads that extend past their mate, overlap, and have clipping on the 3-prime side of one and the 5-prime side of another"
        #[case::overlap_with_mixed_fixed(200, (140, 90), (339, 289), past_mate_clip([25, 0, 0, 175], true), (165, 264, 265, 289))]
        fn clip_pair_past_mate_with_fixed_and_overlap_clipping(
            #[case] read_length: usize,
            #[case] starts: (usize, usize),
            #[case] prior_ends: (usize, usize),
            #[case] clip: Clip,
            #[case] expected: (usize, usize, usize, usize),
        ) {
            let (mut r1, mut r2) = fgbio_pair(read_length, |p| p.start1(starts.0).start2(starts.1));
            assert_eq!((end(&r1), end(&r2)), prior_ends);
            clip_pair_hard(&clip, &mut r1, &mut r2);
            assert_eq!((start(&r1), end(&r1), start(&r2), end(&r2)), expected);
        }

        /// `ClipBamTest.scala:621` "unmap reads when the hard clipping length requested is
        /// greater than the length of the reads". fgbio's `UnmappedStart` (0) corresponds to BAM
        /// `POS` -1 with no alignment start or end.
        #[test]
        fn clip_pair_unmaps_reads_when_clipping_exceeds_read_length() {
            let (mut r1, mut r2) = fgbio_pair(100, |p| p.start1(100).start2(300));
            assert_eq!((end(&r1), end(&r2)), (199, 399));
            clip_pair_hard(&past_mate_clip([101, 0, 101, 0], true), &mut r1, &mut r2);
            assert_eq!((r1.is_unmapped(), r2.is_unmapped()), (true, true));
            assert_eq!((r1.pos(), r2.pos()), (-1, -1));
            assert_eq!((r1.alignment_start_1based(), r2.alignment_start_1based()), (None, None));
            assert_eq!((r1.alignment_end_1based(), r2.alignment_end_1based()), (None, None));
        }
    }
}

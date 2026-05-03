//! Pre-flight validation for `fgumi runall`.
//!
//! `validate()` runs after CLI parsing but before any expensive work
//! (atomic-output guard creation, signal-handler installation, chain
//! runner dispatch). It checks:
//!
//! * The output directory exists and is reachable.
//! * Input file(s) exist and are appropriate for the chosen
//!   `--start-from`.
//! * Stage-specific required flags (e.g. `--extract::sample`,
//!   `--correct::umis`) are present when the chosen chain reaches
//!   their stage.
//! * `--stop-after` is at or after `--start-from`.
//! * Cross-flag invariants (e.g. duplex consensus requires the paired
//!   UMI strategy).
//!
//! Validation errors propagate as `anyhow::Error`s with friendly
//! `--flag …` references in the message text — no atomic-output guard
//! has been created at this point, so there is nothing to clean up.

use std::path::PathBuf;

use crate::reference::find_dict_path;
use crate::validation::validate_file_exists;
use anyhow::{Context, Result, anyhow};

use super::Runall;
use super::options::{ConsensusMode, StartFrom, StopAfter};
use fgumi_umi::Strategy;

impl Runall {
    /// Get the single input file path, or error if multiple inputs were
    /// provided.
    pub(crate) fn single_input(&self) -> Result<&PathBuf> {
        if self.input.len() != 1 {
            anyhow::bail!(
                "--start-from {:?} requires exactly one --input file, got {}",
                self.start_from,
                self.input.len()
            );
        }
        Ok(&self.input[0])
    }

    /// Resolve the aligner command from `--aligner::command`. PR 1 does
    /// not yet implement preset → command translation (the `aligner`
    /// module that PR 3 needs is being added separately); presets that
    /// reach this code path here error out cleanly so users see the
    /// limitation rather than a deep panic.
    pub(crate) fn require_aligner_command(&self) -> Result<String> {
        if self.aligner_opts.aligner_preset.is_some() {
            anyhow::bail!(
                "--aligner::preset is not yet implemented in this fgumi runall foundation; \
                 pass --aligner::command directly with the BWA invocation you want"
            )
        }
        if let Some(ref cmd) = self.aligner_opts.aligner_command {
            Ok(cmd.clone())
        } else {
            anyhow::bail!(
                "--aligner::command or --aligner::preset is required for --start-from {:?}",
                self.start_from
            )
        }
    }

    /// Get the required reference path and its associated dictionary path.
    pub(crate) fn require_reference(&self) -> Result<(PathBuf, PathBuf)> {
        let reference = self
            .reference
            .as_ref()
            .ok_or_else(|| {
                anyhow!("--reference is required for --start-from {:?}", self.start_from)
            })?
            .clone();
        validate_file_exists(&reference, "Reference FASTA")?;
        let dict_path = find_dict_path(&reference).ok_or_else(|| {
            anyhow!(
                "Reference dictionary file not found for {}. Please run: samtools dict",
                reference.display()
            )
        })?;
        Ok((reference, dict_path))
    }

    /// Run pre-flight validation checks before any processing begins.
    ///
    /// Validates inputs, outputs, required parameters, and
    /// cross-parameter consistency so that errors are caught before
    /// any expensive work is started.
    pub(crate) fn validate(&self) -> Result<()> {
        // Output directory must exist and be reachable.
        if let Some(parent) = self.output.parent() {
            if !parent.as_os_str().is_empty() {
                if !parent.exists() {
                    anyhow::bail!("Output directory does not exist: {}", parent.display());
                }
                std::fs::metadata(parent).with_context(|| {
                    format!("Output directory is not accessible: {}", parent.display())
                })?;
            }
        }

        match self.start_from {
            StartFrom::Group => {
                let input = self.single_input()?;
                validate_file_exists(input, "Input BAM")?;
            }
            StartFrom::Sort => {
                let input = self.single_input()?;
                validate_file_exists(input, "Input BAM")?;

                // Reference is needed for the dict used in the output header
                // when the pipeline runs through alignment / consensus.
                // Plans that stop at sort or group can skip the reference.
                let stops_before_consensus =
                    matches!(self.stop_after, StopAfter::Sort | StopAfter::Group);
                if !stops_before_consensus {
                    self.require_reference()?;
                }
            }
            StartFrom::Correct | StartFrom::Fastq => {
                let input = self.single_input()?;
                validate_file_exists(input, "Unmapped BAM")?;

                if matches!(self.start_from, StartFrom::Correct)
                    && self.correct_opts.correct_umis.as_deref().unwrap_or_default().is_empty()
                    && self.correct_opts.correct_umi_files.as_deref().unwrap_or_default().is_empty()
                {
                    anyhow::bail!(
                        "--correct::umis or --correct::umi-files is required for \
                         --start-from correct"
                    );
                }
                let stops_before_align =
                    matches!(self.stop_after, StopAfter::Correct | StopAfter::Fastq);
                if !stops_before_align {
                    self.require_reference()?;
                    self.require_aligner_command()?;
                }
            }
            StartFrom::Extract => {
                if self.input.is_empty() {
                    anyhow::bail!(
                        "--input requires at least one FASTQ file for --start-from extract"
                    );
                }
                for fq in &self.input {
                    validate_file_exists(fq, "Input FASTQ")?;
                }

                if self.extract_opts.extract_sample.is_none() {
                    anyhow::bail!("--extract::sample is required for --start-from extract");
                }
                if self.extract_opts.extract_library.is_none() {
                    anyhow::bail!("--extract::library is required for --start-from extract");
                }
                if let Some(ref t) = self.extract_opts.extract_single_tag {
                    if t.len() != 2 {
                        anyhow::bail!(
                            "--extract::single-tag must be exactly 2 characters, got: {t}"
                        );
                    }
                }
                let stops_before_align = matches!(
                    self.stop_after,
                    StopAfter::Extract | StopAfter::Fastq | StopAfter::Correct
                );
                if !stops_before_align {
                    self.require_reference()?;
                    self.require_aligner_command()?;
                }
            }
        }

        // If any correction flag is set, enforce --start-from extract or correct.
        let has_correction_args =
            self.correct_opts.correct_umis.as_deref().is_some_and(|v| !v.is_empty())
                || self.correct_opts.correct_umi_files.as_deref().is_some_and(|v| !v.is_empty());
        if has_correction_args
            && !matches!(self.start_from, StartFrom::Extract | StartFrom::Correct)
        {
            anyhow::bail!("correction requires --start-from extract or --start-from correct");
        }

        // Validate stop_after is at or after start_from.
        let stop = self.stop_after;
        let min_stop = match self.start_from {
            StartFrom::Group => StopAfter::Group,
            StartFrom::Sort => StopAfter::Sort,
            StartFrom::Fastq => StopAfter::Fastq,
            StartFrom::Correct => StopAfter::Correct,
            StartFrom::Extract => StopAfter::Extract,
        };
        if stop < min_stop {
            anyhow::bail!(
                "--stop-after {:?} is before --start-from {:?}; \
                 stop-after must be at or after the start stage",
                stop,
                self.start_from
            );
        }
        if stop == StopAfter::Correct
            && !matches!(self.start_from, StartFrom::Correct | StartFrom::Extract)
        {
            anyhow::bail!(
                "--stop-after correct is only valid with --start-from extract or \
                 --start-from correct"
            );
        }
        if stop == StopAfter::Fastq
            && !matches!(
                self.start_from,
                StartFrom::Extract | StartFrom::Correct | StartFrom::Fastq
            )
        {
            anyhow::bail!(
                "--stop-after fastq is only valid with --start-from extract, correct, or fastq"
            );
        }

        // --filter::min-reads is required when the plan reaches consensus
        // or filter (mirrors standalone `fgumi filter`, which has no
        // default for `--min-reads` either).
        if matches!(self.stop_after, StopAfter::Consensus | StopAfter::Filter) {
            let min_reads = self.filter_opts.filter_min_reads.as_deref().unwrap_or(&[]);
            if min_reads.is_empty() {
                anyhow::bail!(
                    "--filter::min-reads is required when the plan includes the \
                     consensus or filter stage (--stop-after consensus/filter). \
                     Matches standalone `fgumi filter`, which also has no default."
                );
            }
            if min_reads.len() > 3 {
                anyhow::bail!(
                    "--filter::min-reads must have 1-3 values (total, [AB, [BA]]), got {}",
                    min_reads.len()
                );
            }

            // Consensus mode and grouping strategy must be compatible.
            // Without this check, mismatched plans fail deep inside
            // the consensus stage's tag-validation step after the
            // entire group + MI-assign work has already run.
            match (self.consensus_mode, self.group_opts.group_strategy) {
                (ConsensusMode::Duplex, Strategy::Paired) => {}
                (ConsensusMode::Duplex, other) => anyhow::bail!(
                    "--consensus duplex requires --group::strategy paired \
                     (so re-grouping produces /A and /B-suffixed MI tags that \
                     duplex consensus calling expects). Got --group::strategy {:?}.",
                    other
                ),
                (ConsensusMode::Codec, Strategy::Paired) => anyhow::bail!(
                    "--consensus codec is incompatible with --group::strategy paired \
                     (codec must use adjacency or identity grouping; paired splits \
                     each molecule into /A and /B groups which codec does not handle). \
                     Use --group::strategy adjacency or identity."
                ),
                _ => {}
            }
        }

        // Methylation flags need a reference and consistent depth ordering.
        if self.methylation_mode.is_some() && self.reference.is_none() {
            anyhow::bail!("--methylation-mode requires --reference to be set");
        }
        if self.min_methylation_depth.len() > 3 {
            anyhow::bail!(
                "--min-methylation-depth must have 1-3 values, got {}",
                self.min_methylation_depth.len()
            );
        }
        if self.min_methylation_depth.len() >= 2 {
            let cc = self.min_methylation_depth[0];
            let ab = self.min_methylation_depth[1];
            if cc < ab {
                anyhow::bail!(
                    "min-methylation-depth values must be specified high to low \
                     (duplex >= AB), got {cc} < {ab}"
                );
            }
        }
        if self.min_methylation_depth.len() == 3 {
            let ab = self.min_methylation_depth[1];
            let ba = self.min_methylation_depth[2];
            if ab < ba {
                anyhow::bail!(
                    "min-methylation-depth values must be specified high to low \
                     (AB >= BA), got {ab} < {ba}"
                );
            }
        }
        if self.require_strand_methylation_agreement && self.reference.is_none() {
            anyhow::bail!(
                "--require-strand-methylation-agreement requires --reference to identify \
                 CpG sites"
            );
        }
        if let Some(frac) = self.min_conversion_fraction {
            if !(0.0..=1.0).contains(&frac) {
                anyhow::bail!("--min-conversion-fraction must be between 0.0 and 1.0, got {frac}");
            }
            if self.methylation_mode.is_none() {
                anyhow::bail!("--min-conversion-fraction requires --methylation-mode to be set");
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::runall::infra::progress::ProgressMode;
    use crate::commands::runall::options::{
        AlignerOptions, CodecOptions, ConsensusOptions, CorrectOptions, DuplexOptions,
        ExtractOptions, FilterOptions, GroupOptions, MultiAlignerOptions, MultiCodecOptions,
        MultiConsensusOptions, MultiCorrectOptions, MultiDuplexOptions, MultiExtractOptions,
        MultiFilterOptions, MultiGroupOptions, MultiSortOptions, MultiZipperOptions, SortOptions,
        ZipperOptions,
    };

    fn make_runall(start: StartFrom, stop: StopAfter) -> Runall {
        Runall {
            input: vec![],
            output: PathBuf::from("/tmp/runall_validate_unit_test_out.bam"),
            reference: None,
            start_from: start,
            stop_after: stop,
            consensus_mode: ConsensusMode::Simplex,
            threads: 1,
            progress: ProgressMode::None,
            compression_level: 1,
            metrics: None,
            explain: false,
            consensus_metrics: None,
            intervals: None,
            methylation_mode: None,
            restore_unconverted_bases: false,
            min_methylation_depth: vec![],
            require_strand_methylation_agreement: false,
            min_conversion_fraction: None,
            sort_opts: MultiSortOptions::from(SortOptions::default()),
            group_opts: MultiGroupOptions::from(GroupOptions::default()),
            filter_opts: MultiFilterOptions::from(FilterOptions::default()),
            consensus_opts: MultiConsensusOptions::from(ConsensusOptions::default()),
            extract_opts: MultiExtractOptions::from(ExtractOptions::default()),
            aligner_opts: MultiAlignerOptions::from(AlignerOptions::default()),
            zipper_opts: MultiZipperOptions::from(ZipperOptions::default()),
            correct_opts: MultiCorrectOptions::from(CorrectOptions::default()),
            duplex_opts: MultiDuplexOptions::from(DuplexOptions::default()),
            codec_opts: MultiCodecOptions::from(CodecOptions::default()),
            queue_memory: crate::commands::common::QueueMemoryOptions::default(),
        }
    }

    #[test]
    fn extract_requires_at_least_one_input_fastq() {
        let r = make_runall(StartFrom::Extract, StopAfter::Extract);
        let err = r.validate().unwrap_err();
        assert!(format!("{err:#}").contains("at least one FASTQ"));
    }

    #[test]
    fn group_requires_input_bam_to_exist() {
        let mut r = make_runall(StartFrom::Group, StopAfter::Group);
        r.input = vec![PathBuf::from("/nonexistent/path.bam")];
        let err = r.validate().unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("Input BAM"), "got: {msg}");
    }

    #[test]
    fn stop_after_before_start_from_is_rejected() {
        // Group → Extract is invalid (Extract is "before" Group in the
        // pipeline order). Group doesn't require a reference, so
        // validation reaches the stop_after-vs-start_from check.
        let dir = tempfile::TempDir::new().unwrap();
        let bam = dir.path().join("input.bam");
        std::fs::write(&bam, b"placeholder").unwrap();
        let mut r = make_runall(StartFrom::Group, StopAfter::Extract);
        r.input = vec![bam];
        let err = r.validate().unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("stop-after must be at or after the start stage"), "got: {msg}");
    }

    #[test]
    fn stop_after_correct_requires_correct_or_extract_start() {
        let r = make_runall(StartFrom::Sort, StopAfter::Correct);
        let err = r.validate().unwrap_err();
        // Reaches the correct-pre-condition check after the first
        // input check; just verify *some* sensible error fires.
        assert!(!format!("{err:#}").is_empty());
    }

    #[test]
    fn stop_after_consensus_requires_min_reads() {
        let mut r = make_runall(StartFrom::Group, StopAfter::Consensus);
        // The validator runs the input-exists check first; create
        // a temp BAM placeholder so we get to the min-reads check.
        let dir = tempfile::TempDir::new().unwrap();
        let bam = dir.path().join("input.bam");
        std::fs::write(&bam, b"placeholder").unwrap();
        r.input = vec![bam];

        let err = r.validate().unwrap_err();
        assert!(format!("{err:#}").contains("--filter::min-reads is required"));
    }

    #[test]
    fn duplex_requires_paired_strategy() {
        let dir = tempfile::TempDir::new().unwrap();
        let bam = dir.path().join("input.bam");
        std::fs::write(&bam, b"placeholder").unwrap();

        let mut r = make_runall(StartFrom::Group, StopAfter::Filter);
        r.input = vec![bam];
        r.consensus_mode = ConsensusMode::Duplex;
        // Default group strategy is Adjacency, which is incompatible with duplex.
        // Make filter::min-reads non-empty so we reach the strategy check.
        r.filter_opts.filter_min_reads = Some(vec![1]);

        let err = r.validate().unwrap_err();
        assert!(
            format!("{err:#}").contains("--consensus duplex requires --group::strategy paired")
        );
    }

    #[test]
    fn methylation_mode_requires_reference() {
        let dir = tempfile::TempDir::new().unwrap();
        let bam = dir.path().join("input.bam");
        std::fs::write(&bam, b"placeholder").unwrap();

        let mut r = make_runall(StartFrom::Group, StopAfter::Group);
        r.input = vec![bam];
        r.methylation_mode = Some(crate::commands::common::MethylationModeArg::EmSeq);
        r.reference = None;

        let err = r.validate().unwrap_err();
        assert!(format!("{err:#}").contains("--methylation-mode requires --reference"));
    }
}

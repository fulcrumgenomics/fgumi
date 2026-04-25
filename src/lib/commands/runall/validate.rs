//! Validation logic for the pipeline command.

use std::path::PathBuf;

use crate::reference::find_dict_path;
use crate::validation::validate_file_exists;
use anyhow::{Context, Result, anyhow};

use super::Runall;
use super::options::{StartFrom, StopAfter};

impl Runall {
    /// Get the single input file path, or error if multiple inputs were provided.
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

    /// Resolve the aligner command from `--aligner::command` or `--aligner::preset`.
    ///
    /// If `--aligner::preset` is given, builds the command via [`AlignerPreset::build_command`]
    /// and validates the binary and index files via [`AlignerPreset::validate`]. The reference
    /// must be present (validated separately via [`require_reference`]).
    ///
    /// The returned string may still contain `{ref}` and `{threads}` placeholders if
    /// `--aligner::command` was used directly; callers are responsible for substitution.
    pub(crate) fn require_aligner_command(&self) -> Result<String> {
        if let Some(ref preset_arg) = self.aligner_opts.aligner_preset {
            let reference = self
                .reference
                .as_ref()
                .ok_or_else(|| anyhow!("--reference is required when using --aligner::preset"))?;
            let preset = preset_arg.to_preset();
            preset.validate(reference)?;
            Ok(preset.build_command(
                reference,
                self.aligner_opts.aligner_threads,
                self.aligner_opts.aligner_chunk_size,
                "",
            ))
        } else if let Some(ref cmd) = self.aligner_opts.aligner_command {
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
    /// Validates inputs, outputs, required parameters, and cross-parameter consistency
    /// so that errors are caught before any expensive work is started.
    pub(crate) fn validate(&self) -> Result<()> {
        // Validate output directory exists and is writable
        if let Some(parent) = self.output.parent() {
            // Empty parent means current directory, which always exists
            if !parent.as_os_str().is_empty() && !parent.exists() {
                anyhow::bail!("Output directory does not exist: {}", parent.display());
            }
            if !parent.as_os_str().is_empty() {
                // Check writability by checking metadata (directory must be accessible)
                std::fs::metadata(parent).with_context(|| {
                    format!("Output directory is not accessible: {}", parent.display())
                })?;
            }
        }

        match self.start_from {
            StartFrom::Group => {
                // Requires exactly one input BAM
                let input = self.single_input()?;
                validate_file_exists(input, "Input BAM")?;
            }
            StartFrom::Sort => {
                // Requires exactly one input BAM; reference/dict checked via require_reference
                let input = self.single_input()?;
                validate_file_exists(input, "Input BAM")?;

                // Reference is needed for the dict used in output header when the
                // pipeline runs through the alignment/consensus stages. Plans that
                // stop at sort or group write the input header unchanged and don't
                // need the reference dict.
                let stops_before_consensus =
                    matches!(self.stop_after, StopAfter::Sort | StopAfter::Group);
                if !stops_before_consensus {
                    self.require_reference()?;
                }
            }
            StartFrom::Correct | StartFrom::Fastq => {
                // Requires exactly one unmapped BAM input
                let input = self.single_input()?;
                validate_file_exists(input, "Unmapped BAM")?;

                if matches!(self.start_from, StartFrom::Correct)
                    && self.correct_opts.correct_umis.as_deref().unwrap_or_default().is_empty()
                    && self.correct_opts.correct_umi_files.as_deref().unwrap_or_default().is_empty()
                {
                    anyhow::bail!(
                        "--correct::umis or --correct::umi-files is required for --start-from correct"
                    );
                }
                // Aligner and reference are only required if the pipeline continues past
                // the fastq stage (i.e. not --stop-after correct or --stop-after fastq).
                let stops_before_align =
                    matches!(self.stop_after, StopAfter::Correct | StopAfter::Fastq);
                if !stops_before_align {
                    self.require_reference()?;
                    self.require_aligner_command()?;
                }
            }
            StartFrom::Extract => {
                // Requires FASTQ inputs, sample/library, aligner command, and reference
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
                // Aligner and reference are only required if the pipeline continues past
                // the correct stage.
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

        // Validate stop_after is at or after start_from
        {
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

            // Validate that --stop-after correct is only valid when the flow includes
            // a correction step (i.e. --start-from extract or --start-from correct).
            if stop == StopAfter::Correct
                && !matches!(self.start_from, StartFrom::Correct | StartFrom::Extract)
            {
                anyhow::bail!(
                    "--stop-after correct is only valid with --start-from extract or \
                     --start-from correct"
                );
            }

            // Validate that --stop-after fastq requires a flow that produces FASTQ.
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
        }

        // --filter::min-reads has no default (matches standalone `fgumi filter`,
        // where `--min-reads` is the key parameter the user MUST supply). It is
        // required whenever the plan runs the consensus or filter stage — both
        // depend on it. Plans that stop before consensus (extract / correct /
        // fastq / zipper / sort / group) do not need it.
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
        }

        // Validate methylation options
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
                    "min-methylation-depth values must be specified high to low (duplex >= AB), \
                     got {cc} < {ab}"
                );
            }
        }
        if self.min_methylation_depth.len() == 3 {
            let ab = self.min_methylation_depth[1];
            let ba = self.min_methylation_depth[2];
            if ab < ba {
                anyhow::bail!(
                    "min-methylation-depth values must be specified high to low (AB >= BA), \
                     got {ab} < {ba}"
                );
            }
        }
        if self.require_strand_methylation_agreement && self.reference.is_none() {
            anyhow::bail!(
                "--require-strand-methylation-agreement requires --reference to identify CpG sites"
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

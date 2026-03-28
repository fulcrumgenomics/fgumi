//! Validation logic for the pipeline command.

use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
use fgumi_lib::reference::find_dict_path;
use fgumi_lib::validation::validate_file_exists;

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
            Ok(preset.build_command(reference, self.aligner_opts.aligner_threads, 150_000_000, ""))
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
        // Validate SAM tag lengths (applies to all modes)
        if self.consensus_opts.consensus_mi_tag.len() != 2 {
            anyhow::bail!(
                "MI tag must be exactly 2 characters, got: {}",
                self.consensus_opts.consensus_mi_tag
            );
        }

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
                if self.group_opts.group_umi_tag.len() != 2 {
                    anyhow::bail!(
                        "UMI tag must be exactly 2 characters, got: {}",
                        self.group_opts.group_umi_tag
                    );
                }
                // Reference is needed for the dict used in output header
                self.require_reference()?;
            }
            StartFrom::Zipper => {
                // Requires exactly one mapped input, an unmapped BAM, and reference
                let input = self.single_input()?;
                validate_file_exists(input, "Mapped SAM/BAM")?;
                if self.group_opts.group_umi_tag.len() != 2 {
                    anyhow::bail!(
                        "UMI tag must be exactly 2 characters, got: {}",
                        self.group_opts.group_umi_tag
                    );
                }
                let unmapped_path = self
                    .unmapped
                    .as_ref()
                    .ok_or_else(|| anyhow!("--unmapped is required for --start-from zipper"))?;
                validate_file_exists(unmapped_path, "Unmapped BAM")?;
                self.require_reference()?;
            }
            StartFrom::Correct | StartFrom::Fastq => {
                // Requires exactly one unmapped BAM input
                let input = self.single_input()?;
                validate_file_exists(input, "Unmapped BAM")?;
                if self.group_opts.group_umi_tag.len() != 2 {
                    anyhow::bail!(
                        "UMI tag must be exactly 2 characters, got: {}",
                        self.group_opts.group_umi_tag
                    );
                }
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
                    matches!(self.stop_after, Some(StopAfter::Correct | StopAfter::Fastq));
                if !stops_before_align {
                    self.require_reference()?;
                    self.require_aligner_command()?;
                }
            }
            StartFrom::Align => {
                // Requires FASTQ input, unmapped BAM (for zipper), aligner command or preset,
                // and reference
                let input = self.single_input()?;
                validate_file_exists(input, "Input FASTQ")?;
                if self.group_opts.group_umi_tag.len() != 2 {
                    anyhow::bail!(
                        "UMI tag must be exactly 2 characters, got: {}",
                        self.group_opts.group_umi_tag
                    );
                }
                let unmapped_path = self
                    .unmapped
                    .as_ref()
                    .ok_or_else(|| anyhow!("--unmapped is required for --start-from align"))?;
                validate_file_exists(unmapped_path, "Unmapped BAM")?;
                self.require_reference()?;
                self.require_aligner_command()?;
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
                if self.group_opts.group_umi_tag.len() != 2 {
                    anyhow::bail!(
                        "UMI tag must be exactly 2 characters, got: {}",
                        self.group_opts.group_umi_tag
                    );
                }
                if self.extract_opts.extract_sample.is_none() {
                    anyhow::bail!("--extract::sample is required for --start-from extract");
                }
                if self.extract_opts.extract_library.is_none() {
                    anyhow::bail!("--extract::library is required for --start-from extract");
                }
                if self.extract_opts.extract_cell_tag.len() != 2 {
                    anyhow::bail!(
                        "--extract::cell-tag must be exactly 2 characters, got: {}",
                        self.extract_opts.extract_cell_tag
                    );
                }
                for (name, tag) in [
                    ("--extract::umi-qual-tag", &self.extract_opts.extract_umi_qual_tag),
                    ("--extract::cell-qual-tag", &self.extract_opts.extract_cell_qual_tag),
                    ("--extract::single-tag", &self.extract_opts.extract_single_tag),
                ] {
                    if let Some(t) = tag {
                        if t.len() != 2 {
                            anyhow::bail!("{name} must be exactly 2 characters, got: {t}");
                        }
                    }
                }
                // Aligner and reference are only required if the pipeline continues past
                // the fastq stage.
                let stops_before_align =
                    matches!(self.stop_after, Some(StopAfter::Extract | StopAfter::Fastq));
                if !stops_before_align {
                    self.require_reference()?;
                    self.require_aligner_command()?;
                }
            }
        }

        // Validate stop_after is at or after start_from
        if let Some(stop) = self.stop_after {
            let min_stop = match self.start_from {
                StartFrom::Group => StopAfter::Group,
                StartFrom::Sort => StopAfter::Sort,
                StartFrom::Zipper => StopAfter::Zipper,
                StartFrom::Align => StopAfter::Align,
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
            // a correction step (i.e. --start-from correct).
            if stop == StopAfter::Correct && !matches!(self.start_from, StartFrom::Correct) {
                anyhow::bail!("--stop-after correct is only valid with --start-from correct");
            }

            // Validate that --stop-after fastq requires a flow that produces FASTQ.
            // Currently the FASTQ step is always piped to the aligner, so there's no
            // standalone FASTQ output. Restrict to start-from extract which has the
            // fgumi extract -> fgumi fastq flow.
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

        Ok(())
    }
}

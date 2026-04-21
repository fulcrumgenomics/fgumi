//! Plan construction from `Runall` CLI arguments.
//!
//! [`Runall::build_plan`] is the single dispatch point: given `--start-from`
//! and `--stop-after`, it produces a complete [`Plan`] describing the
//! pipeline plus any post-run finalization work (e.g., correction metrics
//! TSV). Callers:
//!
//! 1. Call [`Runall::build_plan`] to get a [`PlanWithFinalize`].
//! 2. Call [`PlanWithFinalize::run`] to execute end-to-end.
//!
//! All knowledge of "which stages belong where" lives here. The runner
//! ([`super::runner`]) knows only how to instantiate one [`StageSpec`] at a
//! time; it contains no per-entry-point branching.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail, ensure};
use fgumi_umi::TagInfo;
use noodles::sam::Header;

use crate::commands::filter::CollectedFilterMetrics;
use crate::commands::runall::Runall;
use crate::commands::runall::options::{StartFrom, StopAfter};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::driver::PipelineConfig;
use crate::runall::engine::source::fastq::{FastqFormat, detect_fastq_format};
use crate::runall::engine::stages::CorrectMetricsAccumulator;

use super::plan::{
    CorrectionPlanState, ExtractParamsPlanData, Plan, SinkSpec, SourceSpec, StageSpec,
};

// ---------------------------------------------------------------------------
// PlanWithFinalize
// ---------------------------------------------------------------------------

/// A [`Plan`] paired with any out-of-band finalization work that must run
/// after the pipeline completes: UMI-correction metrics TSV, filter stats TSV.
pub(crate) struct PlanWithFinalize {
    pub plan: Plan,
    pub correction_state: Option<CorrectionState>,
    pub filter_stats_state: Option<FilterStatsState>,
}

impl PlanWithFinalize {
    /// Run the plan and finalize metrics/rejects.
    ///
    /// # Errors
    ///
    /// Returns an error if the pipeline fails or if finalization fails.
    pub fn run(self, cancel: &CancelToken) -> Result<()> {
        self.plan.run(cancel)?;
        if let Some(cs) = self.correction_state {
            cs.finalize()?;
        }
        if let Some(fs) = self.filter_stats_state {
            fs.finalize()?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// FilterStatsState
// ---------------------------------------------------------------------------

/// Holds the accumulator + output path for the optional `--filter::stats`
/// TSV. The accumulator is shared with the pipeline's `FilterStage`; finalize
/// folds slots and writes the TSV in the same format as standalone
/// `fgumi filter --stats`.
pub(crate) struct FilterStatsState {
    pub metrics: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
    pub path: PathBuf,
}

impl FilterStatsState {
    /// Fold per-worker slots and write the TSV.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn finalize(self) -> Result<()> {
        let mut total = CollectedFilterMetrics::default();
        for slot in self.metrics.slots() {
            let taken = std::mem::take(&mut *slot.lock());
            total.merge(taken);
        }
        crate::commands::filter::write_filter_stats(
            &self.path,
            total.total_records,
            total.passed_records,
            total.failed_records,
        )?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// CorrectionState (finalization-side mirror of CorrectionPlanState)
// ---------------------------------------------------------------------------

/// Holds all state needed to wire and finalize an optional correction step.
///
/// `CorrectionPlanState` (in [`super::plan`]) is the subset the planner
/// hands to the runner via [`StageSpec::Correct`]. This struct carries the
/// extra fields (umi strings, file paths) needed to write metrics after the
/// pipeline finishes.
pub(crate) struct CorrectionState {
    pub umi_strings: Vec<String>,
    pub umi_length: usize,
    pub encoded: Arc<crate::correct::EncodedUmiSet>,
    pub metrics: Arc<PerThreadAccumulator<CorrectMetricsAccumulator>>,
    pub max_mismatches: usize,
    pub min_distance_diff: usize,
    pub revcomp: bool,
    pub dont_store_original_umis: bool,
    pub cache_size: usize,
    pub rejects_path: Option<PathBuf>,
    pub metrics_path: Option<PathBuf>,
}

impl CorrectionState {
    /// Produce the plan-side subset consumed by [`StageSpec::Correct`].
    pub fn to_plan_state(&self) -> CorrectionPlanState {
        CorrectionPlanState {
            encoded: self.encoded.clone(),
            umi_length: self.umi_length,
            max_mismatches: self.max_mismatches,
            min_distance_diff: self.min_distance_diff,
            revcomp: self.revcomp,
            dont_store_original_umis: self.dont_store_original_umis,
            track_rejects: self.rejects_path.is_some(),
            cache_size: self.cache_size,
            metrics: self.metrics.clone(),
        }
    }

    /// Write metrics (if requested) after the pipeline completes.
    ///
    /// # Errors
    ///
    /// Returns an error if writing the metrics TSV fails.
    pub fn finalize(self) -> Result<()> {
        use fgumi_metrics::correct::UmiCorrectionMetrics;

        // Fold per-worker slots into a single accumulator. The pipeline has
        // returned, so no other writers exist; `into_slots` would suffice,
        // but the per-slot read is cheap and keeps the Arc alive for the
        // length of finalize.
        let mut raw = CorrectMetricsAccumulator::default();
        for slot in self.metrics.slots() {
            let taken = std::mem::take(&mut *slot.lock());
            raw.merge(taken);
        }
        let unmatched_umi = "N".repeat(self.umi_length);

        let mut map: ahash::AHashMap<String, UmiCorrectionMetrics> = ahash::AHashMap::new();
        for umi in self.umi_strings.iter().chain(std::iter::once(&unmatched_umi)) {
            map.entry(umi.clone()).or_insert_with(|| UmiCorrectionMetrics::new(umi.clone()));
        }
        for (umi, counts) in &raw.per_umi {
            let entry =
                map.entry(umi.clone()).or_insert_with(|| UmiCorrectionMetrics::new(umi.clone()));
            entry.total_matches += counts.total_matches;
            entry.perfect_matches += counts.perfect_matches;
            entry.one_mismatch_matches += counts.one_mismatch_matches;
            entry.two_mismatch_matches += counts.two_mismatch_matches;
            entry.other_matches += counts.other_matches;
        }

        let total: u64 = map.values().map(|m| m.total_matches).sum();
        let matched_total: u64 = map
            .iter()
            .filter(|(umi, _)| *umi != &unmatched_umi)
            .map(|(_, m)| m.total_matches)
            .sum();

        #[allow(clippy::cast_precision_loss)]
        for m in map.values_mut() {
            m.fraction_of_matches =
                if total > 0 { m.total_matches as f64 / total as f64 } else { 0.0 };
        }
        let umi_count = map.keys().filter(|u| *u != &unmatched_umi).count();
        #[allow(clippy::cast_precision_loss)]
        let mean = if umi_count > 0 { matched_total as f64 / umi_count as f64 } else { 0.0 };
        for m in map.values_mut() {
            m.representation = if mean > 0.0 { m.total_matches as f64 / mean } else { 0.0 };
        }

        if let Some(path) = self.metrics_path.as_ref() {
            let mut metrics: Vec<UmiCorrectionMetrics> = map.values().cloned().collect();
            metrics.sort_by(|a, b| a.umi.cmp(&b.umi));
            UmiCorrectionMetrics::write_metrics(&metrics, path)?;
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Runall::build_plan — single dispatch point
// ---------------------------------------------------------------------------

impl Runall {
    /// Build a concrete [`Plan`] (plus any finalization work) from the CLI
    /// arguments. Dispatches on `(start_from, stop_after)`.
    ///
    /// # Errors
    ///
    /// Returns an error if required arguments are missing, if input files
    /// cannot be probed (e.g., FASTQ format detection, BAM header read), or
    /// if the requested `(start_from, stop_after)` combination is not yet
    /// wired in `pipeline`.
    pub(crate) fn build_plan(
        &self,
        command_line: &str,
        tmp_output: &Path,
    ) -> Result<PlanWithFinalize> {
        self.build_plan_impl(command_line, tmp_output, false)
    }

    /// Build an explain-mode [`Plan`] without doing any real I/O.
    ///
    /// Unlike [`Self::build_plan`], this method does not open input files,
    /// read BAM headers, load reference dictionaries, or validate aligner
    /// binaries. File-derived fields (BAM headers, reference dict SQ records,
    /// detected FASTQ format, quality encoding) are replaced with stable
    /// placeholders so `--explain` produces a useful plan summary even when
    /// inputs do not exist yet.
    ///
    /// # Errors
    ///
    /// Returns an error only for CLI-level problems that are independent of
    /// file I/O (e.g., missing `--input`, unimplemented `(start_from,
    /// stop_after)` combinations, malformed read structures).
    pub(crate) fn build_plan_explain(
        &self,
        command_line: &str,
        tmp_output: &Path,
    ) -> Result<PlanWithFinalize> {
        self.build_plan_impl(command_line, tmp_output, true)
    }

    fn build_plan_impl(
        &self,
        command_line: &str,
        tmp_output: &Path,
        explain: bool,
    ) -> Result<PlanWithFinalize> {
        match self.start_from {
            StartFrom::Extract => self.plan_from_extract(command_line, tmp_output, explain),
            StartFrom::Correct => self.plan_from_correct(command_line, tmp_output, explain),
            StartFrom::Fastq => self.plan_from_fastq(command_line, tmp_output, explain),
            StartFrom::Sort => self.plan_from_sort(command_line, tmp_output, explain),
            StartFrom::Group => self.plan_from_group(command_line, tmp_output, explain),
        }
    }

    // -----------------------------------------------------------------------
    // start-from: extract
    // -----------------------------------------------------------------------

    fn plan_from_extract(
        &self,
        command_line: &str,
        tmp_output: &Path,
        explain: bool,
    ) -> Result<PlanWithFinalize> {
        ensure!(!self.input.is_empty(), "at least one --input FASTQ is required");

        let fmt = if explain {
            // In explain mode, we cannot probe the input. Default to Plain so the
            // downstream stage prefix choice is deterministic and documented.
            detect_fastq_format(&self.input[0]).unwrap_or(FastqFormat::Plain)
        } else {
            let f = detect_fastq_format(&self.input[0])?;
            for p in &self.input[1..] {
                let other = detect_fastq_format(p)?;
                ensure!(
                    other == f,
                    "all inputs must be the same format; {} is {:?} but {} is {:?}",
                    self.input[0].display(),
                    f,
                    p.display(),
                    other,
                );
            }
            f
        };

        let params: Arc<ExtractParamsPlanData> = Arc::new(self.build_extract_params());
        let read_structures = Arc::new(self.get_read_structures()?);
        let encoding = if explain {
            self.detect_quality_encoding().unwrap_or(crate::extract::QualityEncoding::Standard)
        } else {
            let e = self.detect_quality_encoding()?;
            tracing::info!("Detected quality encoding: {e:?}");
            e
        };

        let unmapped_header = self.build_unmapped_header_or_placeholder(command_line, explain)?;
        let correction_state = self.build_correction_state()?;
        let num_streams = self.input.len();

        let mut stages: Vec<StageSpec> = Vec::new();
        push_fastq_prefix(&mut stages, fmt, num_streams);
        stages.push(StageSpec::Extract { params, read_structures, encoding });
        if let Some(cs) = correction_state.as_ref() {
            stages.push(StageSpec::Correct { state: cs.to_plan_state() });
        }

        let source = SourceSpec::FastqFileRead { paths: self.input.clone() };

        let mut stages = stages;
        let (sink, filter_stats_state, config) = self.append_tail_from_unmapped_bam(
            &mut stages,
            &unmapped_header,
            correction_state.as_ref(),
            command_line,
            tmp_output,
            explain,
            "Pipeline consensus (from extract)",
        )?;

        let (mut stages, mut sink) = (stages, sink);
        wire_deferred_header(&mut stages, &mut sink);
        Ok(PlanWithFinalize {
            plan: Plan { source, stages, sink, config, metrics_tsv: self.metrics.clone() },
            correction_state,
            filter_stats_state,
        })
    }

    // -----------------------------------------------------------------------
    // Shared tail builder: from "unmapped BAM records in the stream" onward.
    //
    // Used by `plan_from_extract` (after Extract + optional Correct),
    // `plan_from_correct` (after Correct), and `plan_from_fastq` (directly from
    // the input BAM). Dispatches on `self.stop_after` and produces the
    // remaining stages, the sink, the filter-stats state, and the per-plan
    // [`PipelineConfig`].
    //
    // `correction_state` is consulted only for StopAfter::Extract | Correct
    // to wire a rejects secondary sink; None otherwise.
    // -----------------------------------------------------------------------

    #[allow(clippy::too_many_arguments)]
    fn append_tail_from_unmapped_bam(
        &self,
        stages: &mut Vec<StageSpec>,
        unmapped_header: &Header,
        correction_state: Option<&CorrectionState>,
        command_line: &str,
        tmp_output: &Path,
        explain: bool,
        consensus_label: &str,
    ) -> Result<(SinkSpec, Option<FilterStatsState>, PipelineConfig)> {
        match self.stop_after {
            StopAfter::Extract | StopAfter::Correct => {
                let (secondary_path, secondary_header) = match correction_state {
                    Some(cs) if cs.rejects_path.is_some() => {
                        (cs.rejects_path.clone(), Some(unmapped_header.clone()))
                    }
                    _ => (None, None),
                };
                append_coalesce_bgzf(stages, self.compression_level);
                Ok((
                    SinkSpec::BamFileWrite {
                        path: tmp_output.to_path_buf(),
                        header: unmapped_header.clone(),
                        secondary_path,
                        secondary_header,
                        queue_capacity: 64,
                    },
                    None,
                    self.default_config(),
                ))
            }
            StopAfter::Fastq => {
                stages.push(StageSpec::ToFastq);
                Ok((
                    SinkSpec::FastqFileSink { path: tmp_output.to_path_buf(), queue_capacity: 64 },
                    None,
                    self.default_config(),
                ))
            }
            StopAfter::Zipper => {
                stages.push(self.align_and_merge_spec_explain(unmapped_header, explain)?);
                append_coalesce_bgzf(stages, self.compression_level);
                let sink_header =
                    self.build_zipper_sink_header_explain(unmapped_header, command_line, explain)?;
                Ok((
                    SinkSpec::BamFileWrite {
                        path: tmp_output.to_path_buf(),
                        header: sink_header,
                        secondary_path: None,
                        secondary_header: None,
                        queue_capacity: 64,
                    },
                    None,
                    self.default_config(),
                ))
            }
            StopAfter::Sort => {
                let sink_header =
                    self.build_zipper_sink_header_explain(unmapped_header, command_line, explain)?;
                stages.push(self.align_and_merge_spec_explain(unmapped_header, explain)?);
                stages.push(self.sort_spec(sink_header.clone()));
                append_coalesce_bgzf(stages, self.compression_level);
                Ok((
                    SinkSpec::BamFileWrite {
                        path: tmp_output.to_path_buf(),
                        header: sink_header,
                        secondary_path: None,
                        secondary_header: None,
                        queue_capacity: 64,
                    },
                    None,
                    self.default_config(),
                ))
            }
            StopAfter::Group => {
                let sink_header =
                    self.build_zipper_sink_header_explain(unmapped_header, command_line, explain)?;
                stages.push(self.align_and_merge_spec_explain(unmapped_header, explain)?);
                stages.push(self.sort_spec(sink_header.clone()));
                stages.push(StageSpec::PositionBatch { header: sink_header.clone() });
                stages.push(self.group_assign_spec());
                stages.push(StageSpec::CoalesceMiGroup {
                    target_chunk_size: default_coalesce_chunk_size(),
                });
                stages.push(StageSpec::BgzfCompress { level: self.compression_level });
                Ok((
                    SinkSpec::BamFileWrite {
                        path: tmp_output.to_path_buf(),
                        header: sink_header,
                        secondary_path: None,
                        secondary_header: None,
                        queue_capacity: 64,
                    },
                    None,
                    self.default_config(),
                ))
            }
            StopAfter::Consensus | StopAfter::Filter => {
                let zipper_header =
                    self.build_zipper_sink_header_explain(unmapped_header, command_line, explain)?;
                let output_header =
                    self.consensus_output_header(&zipper_header, consensus_label, command_line)?;
                let consensus_factory = self.build_consensus_factory(&zipper_header)?;

                stages.push(self.align_and_merge_spec_explain(unmapped_header, explain)?);
                stages.push(self.sort_spec(zipper_header.clone()));
                stages.push(StageSpec::PositionBatch { header: zipper_header.clone() });
                stages.push(self.group_assign_spec());
                stages.push(StageSpec::Consensus {
                    factory: consensus_factory,
                    call_overlapping_bases: self.consensus_opts.consensus_call_overlapping_bases,
                });
                let filter_stats_state = if self.stop_after == StopAfter::Filter {
                    let (spec, state) = self.filter_stage_spec_with_stats(explain)?;
                    stages.push(spec);
                    state
                } else {
                    None
                };
                Ok((
                    SinkSpec::BamSink {
                        path: tmp_output.to_path_buf(),
                        header: output_header,
                        threads: self.threads.max(1),
                        compression_level: self.compression_level,
                        deferred_header: None,
                    },
                    filter_stats_state,
                    self.consensus_config(self.threads.max(1)),
                ))
            }
            StopAfter::Align => bail!(
                "--stop-after align is not supported: the pipeline does not \
                 expose the raw aligner SAM output"
            ),
        }
    }

    // -----------------------------------------------------------------------
    // start-from: correct
    // -----------------------------------------------------------------------

    fn plan_from_correct(
        &self,
        command_line: &str,
        tmp_output: &Path,
        explain: bool,
    ) -> Result<PlanWithFinalize> {
        let input = self.single_input()?;
        let threads = self.threads.max(1);

        let empty_strs: Vec<String> = Vec::new();
        let empty_paths: Vec<PathBuf> = Vec::new();
        let umis_inline = self.correct_opts.correct_umis.as_deref().unwrap_or(&empty_strs);
        let umi_files = self.correct_opts.correct_umi_files.as_deref().unwrap_or(&empty_paths);
        if !explain {
            ensure!(
                !umis_inline.is_empty() || !umi_files.is_empty(),
                "--correct::umis or --correct::umi-files is required for --start-from correct",
            );
        }

        let (umi_strings, umi_length) = if !umis_inline.is_empty() || !umi_files.is_empty() {
            crate::correct::load_umi_sequences(umis_inline, umi_files)?
        } else {
            // Explain-mode fallback: no UMIs provided — use empty placeholder.
            (Vec::<String>::new(), 0)
        };
        let encoded = Arc::new(crate::correct::EncodedUmiSet::new(&umi_strings));
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(self.threads.max(1));
        let cs = CorrectionState {
            umi_strings,
            umi_length,
            encoded,
            metrics,
            max_mismatches: self.correct_opts.correct_max_mismatches,
            min_distance_diff: self.correct_opts.correct_min_distance,
            revcomp: self.correct_opts.correct_revcomp,
            dont_store_original_umis: self.correct_opts.correct_dont_store_original_umis,
            cache_size: self.correct_opts.correct_cache_size,
            rejects_path: self.correct_opts.correct_rejects.clone(),
            metrics_path: self.correct_opts.correct_metrics.clone(),
        };

        let header = if explain {
            let src = crate::runall::engine::source::BamBatchedSource::new(input.clone(), threads);
            src.read_header().unwrap_or_else(|_| Header::default())
        } else {
            let src = crate::runall::engine::source::BamBatchedSource::new(input.clone(), threads);
            src.read_header()?
        };

        let source = SourceSpec::BamBatchedSource { path: input.clone(), threads };
        let mut stages: Vec<StageSpec> = vec![StageSpec::Correct { state: cs.to_plan_state() }];

        let (sink, filter_stats_state, config) = self.append_tail_from_unmapped_bam(
            &mut stages,
            &header,
            Some(&cs),
            command_line,
            tmp_output,
            explain,
            "Pipeline consensus (from correct)",
        )?;

        let (mut stages, mut sink) = (stages, sink);
        wire_deferred_header(&mut stages, &mut sink);
        Ok(PlanWithFinalize {
            plan: Plan { source, stages, sink, config, metrics_tsv: self.metrics.clone() },
            correction_state: Some(cs),
            filter_stats_state,
        })
    }

    // -----------------------------------------------------------------------
    // start-from: fastq
    // -----------------------------------------------------------------------

    fn plan_from_fastq(
        &self,
        command_line: &str,
        tmp_output: &Path,
        explain: bool,
    ) -> Result<PlanWithFinalize> {
        let input = self.single_input()?;
        let threads = self.threads.max(1);

        let header = if explain {
            let src = crate::runall::engine::source::BamBatchedSource::new(input.clone(), threads);
            src.read_header().unwrap_or_else(|_| Header::default())
        } else {
            let src = crate::runall::engine::source::BamBatchedSource::new(input.clone(), threads);
            src.read_header()?
        };

        let source = SourceSpec::BamBatchedSource { path: input.clone(), threads };
        let mut stages: Vec<StageSpec> = Vec::new();

        let (sink, filter_stats_state, config) = self.append_tail_from_unmapped_bam(
            &mut stages,
            &header,
            None,
            command_line,
            tmp_output,
            explain,
            "Pipeline consensus (from fastq)",
        )?;

        let (mut stages, mut sink) = (stages, sink);
        wire_deferred_header(&mut stages, &mut sink);
        Ok(PlanWithFinalize {
            plan: Plan { source, stages, sink, config, metrics_tsv: self.metrics.clone() },
            correction_state: None,
            filter_stats_state,
        })
    }

    // -----------------------------------------------------------------------
    // start-from: sort — runs through consensus + filter to BamSink
    // -----------------------------------------------------------------------

    fn plan_from_sort(
        &self,
        command_line: &str,
        tmp_output: &Path,
        explain: bool,
    ) -> Result<PlanWithFinalize> {
        let input = self.single_input()?;
        let header = if explain {
            crate::bam_io::create_raw_bam_reader(input, 1)
                .map(|(_reader, h)| h)
                .unwrap_or_else(|_| Header::default())
        } else {
            let (_reader, h) = crate::bam_io::create_raw_bam_reader(input, 1)?;
            h
        };
        let output_header =
            self.consensus_output_header(&header, "Pipeline consensus (from sort)", command_line)?;
        let consensus_factory = self.build_consensus_factory(&header)?;
        let threads = self.threads.max(1);

        let source = SourceSpec::BamBatchedSource { path: input.clone(), threads };
        let (filter_spec, filter_stats_state) = self.filter_stage_spec_with_stats(explain)?;
        let stages = vec![
            self.sort_spec(header.clone()),
            StageSpec::PositionBatch { header: header.clone() },
            self.group_assign_spec(),
            StageSpec::Consensus {
                factory: consensus_factory,
                call_overlapping_bases: self.consensus_opts.consensus_call_overlapping_bases,
            },
            filter_spec,
        ];
        let sink = SinkSpec::BamSink {
            path: tmp_output.to_path_buf(),
            header: output_header,
            threads,
            compression_level: self.compression_level,
            deferred_header: None,
        };

        Ok(PlanWithFinalize {
            plan: Plan {
                source,
                stages,
                sink,
                config: self.consensus_config(threads),
                metrics_tsv: self.metrics.clone(),
            },
            correction_state: None,
            filter_stats_state,
        })
    }

    // -----------------------------------------------------------------------
    // start-from: group — same as sort but skips Sort stage
    // -----------------------------------------------------------------------

    fn plan_from_group(
        &self,
        command_line: &str,
        tmp_output: &Path,
        explain: bool,
    ) -> Result<PlanWithFinalize> {
        let input = self.single_input()?;
        let header = if explain {
            crate::bam_io::create_raw_bam_reader(input, 1)
                .map(|(_reader, h)| h)
                .unwrap_or_else(|_| Header::default())
        } else {
            let (_reader, h) = crate::bam_io::create_raw_bam_reader(input, 1)?;
            h
        };
        let output_header =
            self.consensus_output_header(&header, "Pipeline consensus (from group)", command_line)?;
        let consensus_factory = self.build_consensus_factory(&header)?;
        let threads = self.threads.max(1);

        let source = SourceSpec::BamBatchedSource { path: input.clone(), threads };
        let (filter_spec, filter_stats_state) = self.filter_stage_spec_with_stats(explain)?;
        let stages = vec![
            StageSpec::PositionBatch { header: header.clone() },
            self.group_assign_spec(),
            StageSpec::Consensus {
                factory: consensus_factory,
                call_overlapping_bases: self.consensus_opts.consensus_call_overlapping_bases,
            },
            filter_spec,
        ];
        let sink = SinkSpec::BamSink {
            path: tmp_output.to_path_buf(),
            header: output_header,
            threads,
            compression_level: self.compression_level,
            deferred_header: None,
        };

        Ok(PlanWithFinalize {
            plan: Plan {
                source,
                stages,
                sink,
                config: self.consensus_config(threads),
                metrics_tsv: self.metrics.clone(),
            },
            correction_state: None,
            filter_stats_state,
        })
    }

    // -----------------------------------------------------------------------
    // Shared helpers
    // -----------------------------------------------------------------------

    fn default_config(&self) -> PipelineConfig {
        PipelineConfig::auto_tuned(self.threads.max(1))
    }

    /// Tight-queue config used by the consensus path (matches the standalone commands chain).
    fn consensus_config(&self, threads: usize) -> PipelineConfig {
        PipelineConfig {
            worker_threads: threads,
            queue_capacity: 256,
            queue_memory_limit: 128 * 1024 * 1024,
            global_memory_limit: crate::runall::engine::DEFAULT_GLOBAL_MEMORY_LIMIT,
        }
    }

    fn build_correction_state(&self) -> Result<Option<CorrectionState>> {
        let empty_strs: Vec<String> = Vec::new();
        let empty_paths: Vec<PathBuf> = Vec::new();
        let umis_inline = self.correct_opts.correct_umis.as_deref().unwrap_or(&empty_strs);
        let umi_files = self.correct_opts.correct_umi_files.as_deref().unwrap_or(&empty_paths);
        if umis_inline.is_empty() && umi_files.is_empty() {
            return Ok(None);
        }

        let (umi_strings, umi_length) = crate::correct::load_umi_sequences(umis_inline, umi_files)?;
        let encoded = Arc::new(crate::correct::EncodedUmiSet::new(&umi_strings));
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(self.threads.max(1));
        Ok(Some(CorrectionState {
            umi_strings,
            umi_length,
            encoded,
            metrics,
            max_mismatches: self.correct_opts.correct_max_mismatches,
            min_distance_diff: self.correct_opts.correct_min_distance,
            revcomp: self.correct_opts.correct_revcomp,
            dont_store_original_umis: self.correct_opts.correct_dont_store_original_umis,
            cache_size: self.correct_opts.correct_cache_size,
            rejects_path: self.correct_opts.correct_rejects.clone(),
            metrics_path: self.correct_opts.correct_metrics.clone(),
        }))
    }

    /// `align_and_merge_spec` with an explain-mode fallback.
    ///
    /// In explain mode, missing aligner command / reference / dict are
    /// replaced with obvious `(not set)` placeholder strings and paths so
    /// that plan construction can proceed without real files.
    fn align_and_merge_spec_explain(
        &self,
        unmapped_header: &Header,
        explain: bool,
    ) -> Result<StageSpec> {
        let aligner_command = match self.require_aligner_command() {
            Ok(c) => c,
            Err(e) if explain => format!("(not set: {e})"),
            Err(e) => return Err(e),
        };
        let (reference, dict_path) = match self.require_reference() {
            Ok(pair) => pair,
            Err(_) if explain => (PathBuf::from("(not set)"), PathBuf::from("(not set)")),
            Err(e) => return Err(e),
        };
        let tag_info = TagInfo::new(
            self.zipper_opts.zipper_tags_to_remove.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_reverse.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_revcomp.clone().unwrap_or_default(),
        );
        Ok(StageSpec::AlignAndMerge {
            aligner_command,
            reference,
            dict_path,
            unmapped_header: unmapped_header.clone(),
            tag_info,
            skip_pa_tags: self.zipper_opts.zipper_skip_pa_tags,
            no_read_suffix: false,
            deferred_header: None,
        })
    }

    /// [`build_zipper_sink_header`] with an explain-mode fallback. In explain
    /// mode, a missing reference dict yields an empty placeholder header
    /// rather than an error.
    fn build_zipper_sink_header_explain(
        &self,
        unmapped: &Header,
        command_line: &str,
        explain: bool,
    ) -> Result<Header> {
        match self.require_reference() {
            Ok((_reference, dict_path)) => {
                match build_zipper_sink_header(unmapped, &dict_path, command_line) {
                    Ok(h) => Ok(h),
                    Err(_) if explain => Ok(unmapped.clone()),
                    Err(e) => Err(e),
                }
            }
            Err(_) if explain => Ok(unmapped.clone()),
            Err(e) => Err(e),
        }
    }

    /// [`Self::build_unmapped_header`] with an explain-mode placeholder. The
    /// real method requires `--extract::sample` and `--extract::library`; in
    /// explain mode we fall back to a default header if either is unset.
    fn build_unmapped_header_or_placeholder(
        &self,
        command_line: &str,
        explain: bool,
    ) -> Result<Header> {
        if explain
            && (self.extract_opts.extract_sample.is_none()
                || self.extract_opts.extract_library.is_none())
        {
            return Ok(Header::default());
        }
        self.build_unmapped_header(command_line)
    }

    fn sort_spec(&self, header: Header) -> StageSpec {
        StageSpec::Sort {
            memory_limit: self.sort_opts.sort_memory_limit,
            temp_dir: self.sort_opts.sort_temp_dir.clone(),
            threads: self.threads.max(1),
            header,
        }
    }

    fn group_assign_spec(&self) -> StageSpec {
        StageSpec::GroupAssign {
            umi_tag: *crate::sam::SamTag::RX,
            assign_tag: *crate::sam::SamTag::MI,
            strategy: self.group_opts.group_strategy,
            max_edits: self.group_opts.group_max_edits,
            filter_config: self.build_group_filter_config(),
        }
    }

    fn consensus_output_header(
        &self,
        input_header: &Header,
        label: &str,
        command_line: &str,
    ) -> Result<Header> {
        use crate::commands::consensus_runner::create_unmapped_consensus_header;
        create_unmapped_consensus_header(
            input_header,
            &self.consensus_opts.consensus_read_group_id,
            label,
            command_line,
        )
    }
}

// ---------------------------------------------------------------------------
// Module-private helpers
// ---------------------------------------------------------------------------

/// Push the format-dependent FASTQ-decode prefix stages onto `stages`.
///
/// - Bgzf inputs use the raw-block path (`BgzfDecompress` → `FastqBlockParse` → `FastqBlockMerge`).
/// - Gzip/plain inputs use the record-aligned path (`FastqParse` → `FastqPair`).
fn push_fastq_prefix(stages: &mut Vec<StageSpec>, fmt: FastqFormat, num_streams: usize) {
    match fmt {
        FastqFormat::Bgzf => {
            stages.push(StageSpec::BgzfDecompress);
            stages.push(StageSpec::FastqBlockParse);
            stages.push(StageSpec::FastqBlockMerge { num_streams });
        }
        FastqFormat::Gzip | FastqFormat::Plain => {
            stages.push(StageSpec::FastqParse);
            stages.push(StageSpec::FastqPair { num_streams });
        }
    }
}

/// Append the standard `Coalesce<SerializedBatch>` → `BgzfCompress` tail used
/// by every BamFileWrite-terminating pipeline (except the MI-group path,
/// which uses `CoalesceMiGroup` instead).
fn append_coalesce_bgzf(stages: &mut Vec<StageSpec>, compression_level: u32) {
    stages.push(StageSpec::CoalesceSerialized { target_chunk_size: default_coalesce_chunk_size() });
    stages.push(StageSpec::BgzfCompress { level: compression_level });
}

fn default_coalesce_chunk_size() -> usize {
    crate::runall::engine::stages::coalesce::DEFAULT_TARGET_CHUNK_SIZE
}

/// Build the output header for the zipper-merged BAM sink.
///
/// Takes @SQ from the reference dictionary, @RG and @CO from the unmapped
/// header, and appends a runall @PG record.
/// Pair up `StageSpec::AlignAndMerge` with the final `SinkSpec::BamSink`
/// via a shared `Arc<OnceLock<Header>>`, so the sink can use the aligner-
/// aware output header instead of the plan-time approximation built from
/// the dict + unmapped header.
///
/// No-op when the plan does not contain both a compatible stage and sink.
fn wire_deferred_header(stages: &mut [StageSpec], sink: &mut SinkSpec) {
    use std::sync::{Arc, OnceLock};

    // Locate the AlignAndMerge stage's deferred_header slot.
    let Some(stage_slot) = stages.iter_mut().find_map(|s| match s {
        StageSpec::AlignAndMerge { deferred_header, .. } => Some(deferred_header),
        _ => None,
    }) else {
        return;
    };
    // Only wire when the sink is the matching BamSink consensus sink.
    let SinkSpec::BamSink { deferred_header: sink_slot, .. } = sink else {
        return;
    };
    let shared: crate::runall::engine::sink::DeferredHeader = Arc::new(OnceLock::new());
    *stage_slot = Some(shared.clone());
    *sink_slot = Some(shared);
}

fn build_zipper_sink_header(
    unmapped: &Header,
    dict_path: &Path,
    command_line: &str,
) -> Result<Header> {
    let dict_file =
        std::fs::File::open(dict_path).context("Failed to open reference dictionary")?;
    let mut dict_reader = noodles::sam::io::Reader::new(std::io::BufReader::new(dict_file));
    let dict_header = dict_reader.read_header()?;

    let mut builder = Header::builder();
    for (name, map) in dict_header.reference_sequences() {
        builder = builder.add_reference_sequence(name.clone(), map.clone());
    }
    for (id, rg) in unmapped.read_groups() {
        builder = builder.add_read_group(id.clone(), rg.clone());
    }
    for (id, pg) in unmapped.programs().as_ref() {
        builder = builder.add_program(id.clone(), pg.clone());
    }
    for comment in unmapped.comments() {
        builder = builder.add_comment(comment.clone());
    }
    if let Some(hdr) = unmapped.header() {
        builder = builder.set_header(hdr.clone());
    }

    let header = builder.build();
    crate::commands::common::add_pg_record(header, command_line)
}

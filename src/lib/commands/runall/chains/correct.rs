//! `--start-from correct --stop-after correct` chain runner.
//!
//! Thin BAM-linear entry point: opens the input BAM via the same helper
//! the standalone `fgumi correct` command uses, builds
//! [`CorrectPipelineParams`], and forwards to
//! [`crate::commands::correct::run_correct_pipeline`]. The pipeline body
//! (closures, per-thread metrics accumulator, rejects writer,
//! `min_corrected` check, metrics file finalisation) lives exactly once,
//! in `crate::commands::correct`.
//!
//! Gate: matches `(Correct, Correct)` whenever a UMI source is provided
//! (inline UMIs or a UMI file). The `--correct::metrics` and
//! `--correct::min-corrected` flags are forwarded directly to
//! [`run_correct_pipeline`] (which already implements both). The
//! `--correct::rejects` flag remains unsupported pending upstream
//! refactors (fulcrumgenomics/fgumi#329 +
//! fulcrumgenomics/fgumi#330) — when set, the runner still matches the
//! chain shape so the user sees a specific error rather than the generic
//! "not yet implemented" fallthrough. The top-level `--metrics` flag is
//! also honoured: the runner times its `run_correct_pipeline` call and
//! writes a single `MetricsRow { stage: "correct", ... }` row to the
//! supplied TSV path.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use crate::commands::common::{
    CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions, add_pg_record,
};
use crate::commands::correct::{CorrectPipelineParams, EncodedUmiSet, run_correct_pipeline};
use crate::commands::runall::Runall;
use crate::commands::runall::dispatch::{ChainContext, ChainRunner, DispatchContext};
use crate::commands::runall::infra::metrics::{MetricsRow, MetricsWriter};
use crate::commands::runall::options::{StartFrom, StopAfter};
use fgumi_bam_io::create_bam_reader_for_pipeline_with_opts;

/// Stable runner name surfaced via `--explain` and debug logs.
pub(crate) const CORRECT_CHAIN_RUNNER_NAME: &str = "CorrectChainRunner";

/// Chain runner for `--start-from correct --stop-after correct`.
///
/// Built once per `Runall::execute` call; owns the cloned subset of
/// `Runall` fields needed to drive [`run_correct_pipeline`]. Cloning keeps
/// the runner `'static` (so it can be pushed into the trait-object
/// registry) without leaking `Runall` lifetimes.
#[allow(clippy::struct_excessive_bools)] // correction toggles
pub(crate) struct CorrectChainRunner {
    input: PathBuf,
    threads: usize,
    compression_level: u32,
    umis: Vec<String>,
    umi_files: Vec<PathBuf>,
    max_mismatches: usize,
    min_distance_diff: usize,
    cache_size: usize,
    revcomp: bool,
    dont_store_original_umis: bool,
    /// `--correct::rejects`. Still unsupported at the runall layer pending
    /// fulcrumgenomics/fgumi#329 (migrate hand-rolled rejects writers to
    /// first-class `secondary_output`). When set, [`run_chain`] bails with
    /// a specific error pointing at that issue rather than silently
    /// falling through to the not-implemented runner.
    rejects: Option<PathBuf>,
    /// `--correct::metrics`. Forwarded to [`run_correct_pipeline`] via
    /// [`CorrectPipelineParams::metrics_path`]; the per-UMI metrics TSV
    /// is written by [`finalize_correct_metrics`] inside the pipeline.
    correct_metrics: Option<PathBuf>,
    /// `--correct::min-corrected`. Forwarded to [`run_correct_pipeline`]
    /// via [`CorrectPipelineParams::min_corrected`]; the post-pipeline
    /// kept/total-records ratio check is performed there.
    correct_min_corrected: Option<f64>,
    /// `--queue-memory` / `--queue-memory-per-thread` from the top-level
    /// `Runall`. Cloned (rather than per-field-decomposed) so the runner
    /// can call `QueueMemoryOptions::calculate_memory_limit` and
    /// `QueueMemoryOptions::log_memory_config` exactly the way the
    /// standalone `fgumi correct` command does, keeping memory-budget
    /// behaviour in lock-step.
    queue_memory: QueueMemoryOptions,
}

impl CorrectChainRunner {
    /// Build the runner from a parsed [`Runall`] invocation.
    ///
    /// Builds unconditionally: the registry consults `supports()` to decide
    /// whether to dispatch here, so a `Runall` invocation that does not
    /// match `(Correct, Correct)` is allowed to flow through this
    /// constructor. Required-when-running fields are stored as-is and
    /// re-checked inside `run_chain()` so a missing flag surfaces as a
    /// clean `anyhow::Error` rather than a panic.
    pub(crate) fn from_runall(runall: &Runall) -> Self {
        let opts = &runall.correct_opts;
        Self {
            input: runall.input.first().cloned().unwrap_or_default(),
            threads: runall.threads,
            compression_level: runall.compression_level,
            umis: opts.correct_umis.clone().unwrap_or_default(),
            umi_files: opts.correct_umi_files.clone().unwrap_or_default(),
            max_mismatches: opts.correct_max_mismatches,
            min_distance_diff: opts.correct_min_distance,
            cache_size: opts.correct_cache_size,
            revcomp: opts.correct_revcomp,
            dont_store_original_umis: opts.correct_dont_store_original_umis,
            rejects: opts.correct_rejects.clone(),
            correct_metrics: opts.correct_metrics.clone(),
            correct_min_corrected: opts.correct_min_corrected,
            queue_memory: runall.queue_memory.clone(),
        }
    }

    /// Run the chain end-to-end against the prepared
    /// [`ChainContext::tmp_output`] path.
    ///
    /// When `metrics_path` is `Some`, writes a single
    /// [`MetricsRow`] tagged `"correct"` with wall-clock duration and
    /// `records_in == records_out` (the total record count returned by
    /// [`run_correct_pipeline`], which equals kept + rejected). Once
    /// `--correct::rejects` lands (fulcrumgenomics/fgumi#329 +
    /// fulcrumgenomics/fgumi#330), this row should split `records_in` /
    /// `records_out` so the dropped count is visible in the top-level
    /// metrics TSV.
    fn run_chain(
        &self,
        command_line: &str,
        output: &Path,
        metrics_path: Option<&Path>,
    ) -> Result<()> {
        anyhow::ensure!(
            !self.input.as_os_str().is_empty(),
            "CorrectChainRunner: --input is required",
        );
        anyhow::ensure!(
            !self.umis.is_empty() || !self.umi_files.is_empty(),
            "CorrectChainRunner: --correct::umis or --correct::umi-files is required",
        );
        if self.rejects.is_some() {
            // Wiring `--correct::rejects` through the runall correct chain
            // requires the standalone command's hand-rolled
            // `Arc<Mutex<Option<RawBamWriter>>>` to first migrate to the
            // unified pipeline's `secondary_output` (issue #329); the fused
            // extract→correct chain additionally needs FASTQ-pipeline
            // secondary-output support (issue #330's phase 1). Until both
            // land we fail with a specific error so users get a pointer to
            // the tracking work rather than the generic "not yet
            // implemented" fallthrough.
            bail!(
                "--correct::rejects is not yet supported by the runall correct chain; \
                 tracked in fulcrumgenomics/fgumi#329 and fulcrumgenomics/fgumi#330"
            );
        }

        // Load the UMI list from inline strings + files, deduped + uppercased,
        // matching the standalone `fgumi correct` loader.
        let (umi_sequences, umi_length) = load_umi_sequences(&self.umis, &self.umi_files)?;
        let encoded_umi_set = EncodedUmiSet::new(&umi_sequences);

        // Open input BAM with the same streaming-capable helper standalone
        // correct uses, then stamp a fresh `@PG` record on the header.
        let opts = fgumi_bam_io::PipelineReaderOpts::default();
        let (reader, header) = create_bam_reader_for_pipeline_with_opts(&self.input, opts)?;
        let header = add_pg_record(header, command_line)?;

        // Build minimal option structs (defaults match standalone defaults; the
        // few fields the chain runner exposes are passed explicitly via params).
        // `--queue-memory` / `--queue-memory-per-thread` are forwarded from the
        // top-level `Runall` so the chain runner builds the same
        // `BamPipelineConfig.queue_memory_limit` standalone `fgumi correct`
        // would for the same flag values.
        let scheduler_opts = SchedulerOptions::default();
        let compression = CompressionOptions { compression_level: self.compression_level };
        let threading = ThreadingOptions { threads: Some(self.threads) };

        let start = std::time::Instant::now();
        let total_records = run_correct_pipeline(CorrectPipelineParams {
            num_threads: self.threads,
            reader,
            header,
            encoded_umi_set: Arc::new(encoded_umi_set),
            umi_length,
            track_rejects: false,
            scheduler_opts: &scheduler_opts,
            compression: &compression,
            queue_memory: &self.queue_memory,
            threading: &threading,
            rejects_path: None,
            output_path: output,
            metrics_path: self.correct_metrics.as_deref(),
            max_mismatches: self.max_mismatches,
            min_distance_diff: self.min_distance_diff,
            revcomp: self.revcomp,
            cache_size: self.cache_size,
            dont_store_original_umis: self.dont_store_original_umis,
            min_corrected: self.correct_min_corrected,
        })
        .with_context(|| format!("CorrectChainRunner pipeline writing to {}", output.display()))?;
        let wall_time_secs = start.elapsed().as_secs_f64();

        if let Some(path) = metrics_path {
            // Until `--correct::rejects` is wired (issues #329 + #330) we
            // can't separate the kept count from the total here, so both
            // columns carry the total record count returned by
            // `run_correct_pipeline`.
            let mut writer = MetricsWriter::create(path)
                .with_context(|| format!("opening --metrics path {} for write", path.display()))?;
            writer
                .write_row(&MetricsRow {
                    stage: "correct".into(),
                    wall_time_secs,
                    records_in: total_records,
                    records_out: total_records,
                })
                .with_context(|| format!("writing correct metrics row to {}", path.display()))?;
        }

        Ok(())
    }
}

impl ChainRunner for CorrectChainRunner {
    fn name(&self) -> &'static str {
        CORRECT_CHAIN_RUNNER_NAME
    }

    /// Match `(Correct, Correct)` whenever a UMI source is provided
    /// (inline `--correct::umis` or `--correct::umi-files`).
    ///
    /// `--correct::metrics` and `--correct::min-corrected` are forwarded
    /// directly to [`run_correct_pipeline`] (which already owns both
    /// behaviours), so they don't gate the runner.
    ///
    /// `--correct::rejects` ALSO does not gate `supports()` — when set we
    /// still want this runner to win over `NotImplementedRunner` so the
    /// user sees the specific "tracked in #329 / #330" error from
    /// [`run_chain`] instead of the generic fallthrough message.
    ///
    /// The top-level `--metrics` flag does NOT gate this runner — it is
    /// honoured directly by `run_chain`, which writes a single
    /// `correct` row to the supplied TSV.
    fn supports(&self, ctx: &DispatchContext<'_>) -> bool {
        if self.umis.is_empty() && self.umi_files.is_empty() {
            return false;
        }
        ctx.start_from == StartFrom::Correct && ctx.stop_after == StopAfter::Correct
    }

    fn run(&self, ctx: ChainContext<'_>) -> Result<()> {
        let metrics_path = ctx.dispatch.metrics_path;
        self.run_chain(ctx.command_line, &ctx.tmp_output, metrics_path)
    }
}

/// Load and deduplicate UMI sequences from inline strings and/or files.
///
/// Mirrors `CorrectUmis::load_umi_sequences`; reproduced here so the chain
/// runner doesn't need to construct a full `CorrectUmis`. All UMIs must
/// share a single length.
pub(crate) fn load_umi_sequences(
    umis: &[String],
    umi_files: &[std::path::PathBuf],
) -> Result<(Vec<String>, usize)> {
    let mut umi_set: std::collections::HashSet<String> =
        umis.iter().map(|s| s.to_uppercase()).collect();

    for file in umi_files {
        let content = std::fs::read_to_string(file)
            .with_context(|| format!("read UMI file {}", file.display()))?;
        for line in content.lines() {
            let umi = line.trim().to_uppercase();
            if !umi.is_empty() {
                umi_set.insert(umi);
            }
        }
    }

    if umi_set.is_empty() {
        anyhow::bail!("No UMIs provided.");
    }

    let mut umi_sequences: Vec<String> = umi_set.into_iter().collect();
    umi_sequences.sort_unstable();

    let first_len = umi_sequences[0].len();
    if !umi_sequences.iter().all(|u| u.len() == first_len) {
        anyhow::bail!("All UMIs must have the same length.");
    }

    Ok((umi_sequences, first_len))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::runall::infra::progress::ProgressMode;
    use crate::commands::runall::options::{
        AlignerOptions, CodecOptions, ConsensusMode, ConsensusOptions, CorrectOptions,
        DuplexOptions, ExtractOptions, FilterOptions, GroupOptions, MultiAlignerOptions,
        MultiCodecOptions, MultiConsensusOptions, MultiCorrectOptions, MultiDuplexOptions,
        MultiExtractOptions, MultiFilterOptions, MultiGroupOptions, MultiSortOptions,
        MultiZipperOptions, SortOptions, ZipperOptions,
    };
    use std::path::PathBuf;

    fn runall_for_correct(input: PathBuf, output: PathBuf, umis: Vec<String>) -> Runall {
        let correct = CorrectOptions { umis: Some(umis), ..CorrectOptions::default() };
        Runall {
            input: vec![input],
            output,
            reference: None,
            start_from: StartFrom::Correct,
            stop_after: StopAfter::Correct,
            consensus_mode: ConsensusMode::Simplex,
            threads: 4,
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
            correct_opts: MultiCorrectOptions::from(correct),
            duplex_opts: MultiDuplexOptions::from(DuplexOptions::default()),
            codec_opts: MultiCodecOptions::from(CodecOptions::default()),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    fn dispatch_ctx(start: StartFrom, stop: StopAfter) -> DispatchContext<'static> {
        DispatchContext { start_from: start, stop_after: stop, metrics_path: None }
    }

    #[test]
    fn supports_only_correct_to_correct_with_umi_source() {
        let dir = tempfile::TempDir::new().unwrap();
        let runall = runall_for_correct(
            dir.path().join("in.bam"),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        let runner = CorrectChainRunner::from_runall(&runall);

        assert!(runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Correct)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Extract)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Correct)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Filter)));

        let metrics_path = dir.path().join("metrics.tsv");
        let with_metrics = DispatchContext {
            start_from: StartFrom::Correct,
            stop_after: StopAfter::Correct,
            metrics_path: Some(metrics_path.as_path()),
        };
        assert!(
            runner.supports(&with_metrics),
            "correct chain runner now writes its own --metrics row; the registry \
             must dispatch to it even when --metrics is set"
        );
    }

    #[test]
    fn supports_when_correct_rejects_set_so_run_can_emit_specific_error() {
        // Lifted gate (this PR): we want the runner to claim the chain
        // shape even when `--correct::rejects` is set, so the user sees
        // the specific "tracked in #329 / #330" error from `run_chain`
        // rather than the generic NotImplementedRunner fallthrough.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_correct(
            dir.path().join("in.bam"),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.correct_opts.correct_rejects = Some(dir.path().join("rejects.bam"));
        let runner = CorrectChainRunner::from_runall(&runall);
        assert!(
            runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Correct)),
            "--correct::rejects must NOT gate supports() — it routes through run_chain to fail \
             with a specific error pointing at issues #329 and #330"
        );
    }

    #[test]
    fn supports_with_correct_metrics_set() {
        // `--correct::metrics` is now forwarded directly to
        // `run_correct_pipeline`; the runner should claim the chain shape.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_correct(
            dir.path().join("in.bam"),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.correct_opts.correct_metrics = Some(dir.path().join("correct-metrics.tsv"));
        let runner = CorrectChainRunner::from_runall(&runall);
        assert!(
            runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Correct)),
            "--correct::metrics must no longer gate supports() (this PR)"
        );
    }

    #[test]
    fn supports_with_correct_min_corrected_set() {
        // `--correct::min-corrected` is now forwarded directly to
        // `run_correct_pipeline`; the runner should claim the chain shape.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_correct(
            dir.path().join("in.bam"),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.correct_opts.correct_min_corrected = Some(0.5);
        let runner = CorrectChainRunner::from_runall(&runall);
        assert!(
            runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Correct)),
            "--correct::min-corrected must no longer gate supports() (this PR)"
        );
    }

    #[test]
    fn rejects_set_returns_specific_error_pointing_at_issues_329_330() {
        // Lifted gate (this PR) means `supports()` is now true; the
        // specific error must surface from `run_chain` itself.
        let dir = tempfile::TempDir::new().unwrap();
        let in_bam = dir.path().join("in.bam");
        write_test_bam(&in_bam, 1, &[b"AAAAAAAA"]);
        let mut runall = runall_for_correct(
            in_bam.clone(),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.correct_opts.correct_rejects = Some(dir.path().join("rejects.bam"));
        let runner = CorrectChainRunner::from_runall(&runall);
        let err = runner
            .run_chain("fgumi runall", &dir.path().join("ignored.bam"), None)
            .expect_err("run_chain must bail when --correct::rejects is set");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("fulcrumgenomics/fgumi#329"),
            "error must mention upstream issue #329; got: {msg}"
        );
        assert!(
            msg.contains("fulcrumgenomics/fgumi#330"),
            "error must mention upstream issue #330; got: {msg}"
        );
    }

    #[test]
    fn from_runall_propagates_queue_memory_options() {
        // Top-level `--queue-memory` / `--queue-memory-per-thread` must
        // reach the chain runner so it can build the same
        // `BamPipelineConfig.pipeline.queue_memory_limit` standalone
        // `fgumi correct` would compute for the same values.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_correct(
            dir.path().join("in.bam"),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.queue_memory = QueueMemoryOptions {
            queue_memory: "2GB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let runner = CorrectChainRunner::from_runall(&runall);
        assert_eq!(runner.queue_memory.queue_memory, "2GB");
        assert!(!runner.queue_memory.queue_memory_per_thread);
    }

    /// Build a tiny test BAM at `path` containing `n` single-segment
    /// records, each carrying an 8 bp UMI in the RX tag from `umi_cycle`.
    /// Mirrors the `create_test_bam` helper in `commands::correct::tests`
    /// (kept local so this module doesn't depend on its private test
    /// utilities).
    fn write_test_bam(path: &Path, n: usize, umi_cycle: &[&[u8]]) {
        use fgumi_raw_bam::{
            SamBuilder as RawSamBuilder, raw_record_to_record_buf, testutil::encode_op,
        };
        use noodles::sam::Header;
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::header::record::value::map::Map;

        let header = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZero::new(1000).expect("non-zero reference length"),
                ),
            )
            .build();

        let mut writer = noodles::bam::io::writer::Builder.build_from_path(path).unwrap();
        writer.write_header(&header).unwrap();
        for i in 0..n {
            let mut b = RawSamBuilder::new();
            let name = format!("r{i}");
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(0)
                .mapq(60)
                .cigar_ops(&[encode_op(0, 10)])
                .sequence(b"AAAAAAAAAA")
                .qualities(&[40u8; 10]);
            b.add_string_tag(crate::sam::SamTag::RX, umi_cycle[i % umi_cycle.len()]);
            let raw = b.build();
            let record = raw_record_to_record_buf(&raw, &noodles::sam::Header::default()).unwrap();
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish().unwrap();
    }

    #[test]
    fn metrics_tsv_emits_one_correct_row() {
        // Wire `--metrics` end-to-end through `run_chain` with a real BAM
        // input, then parse the TSV. Verifies header, single data row,
        // stable column ordering, and `records_in == records_out` (true
        // for the supported chain configurations until PR #329 lands).
        let dir = tempfile::TempDir::new().unwrap();
        let in_bam = dir.path().join("in.bam");
        let out_tmp = dir.path().join("out.bam.tmp");
        let metrics = dir.path().join("metrics.tsv");

        // 8 records, all UMIs in the allowed list — every record passes
        // through, no rejections.
        let umis: &[&[u8]] = &[b"AAAAAAAA", b"CCCCCCCC", b"GGGGGGGG", b"TTTTTTTT"];
        write_test_bam(&in_bam, 8, umis);

        let runall = runall_for_correct(
            in_bam.clone(),
            dir.path().join("ignored.bam"),
            vec![
                "AAAAAAAA".to_string(),
                "CCCCCCCC".to_string(),
                "GGGGGGGG".to_string(),
                "TTTTTTTT".to_string(),
            ],
        );
        let runner = CorrectChainRunner::from_runall(&runall);
        runner
            .run_chain("fgumi runall", &out_tmp, Some(metrics.as_path()))
            .expect("chain runs with --metrics");

        let contents = std::fs::read_to_string(&metrics).expect("read metrics");
        let mut lines = contents.lines();
        assert_eq!(
            lines.next().unwrap(),
            "stage\twall_time_secs\trecords_in\trecords_out",
            "header row"
        );
        let row = lines.next().expect("one correct row");
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), 4, "row has 4 columns: {row}");
        assert_eq!(cols[0], "correct", "stage column");
        let wall: f64 = cols[1].parse().expect("wall_time_secs is a float");
        assert!(wall >= 0.0, "wall_time_secs is non-negative: {wall}");
        assert_eq!(cols[2], "8", "records_in: 8 input records");
        assert_eq!(
            cols[3], "8",
            "records_out: 8 (records_in == records_out without --correct::rejects/min-corrected)"
        );
        assert!(lines.next().is_none(), "exactly one data row");
    }

    #[test]
    fn build_pipeline_config_applies_queue_memory_limit() {
        // End-to-end check that `--queue-memory` reaches `BamPipelineConfig`:
        // 1 GiB total (per-thread = false) → exactly 1 GiB queue-memory limit
        // on the resulting config, regardless of `--threads`.
        use crate::commands::common::build_pipeline_config;
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_correct(
            dir.path().join("in.bam"),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.threads = 4;
        runall.queue_memory = QueueMemoryOptions {
            queue_memory: "1GiB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let runner = CorrectChainRunner::from_runall(&runall);
        // Reconstruct the same pipeline config `run_chain` builds so we can
        // assert on the resolved limit without spinning up a real BAM.
        let scheduler_opts = SchedulerOptions::default();
        let compression = CompressionOptions { compression_level: runner.compression_level };
        let config = build_pipeline_config(
            &scheduler_opts,
            &compression,
            &runner.queue_memory,
            runner.threads,
        )
        .expect("config builds");
        assert_eq!(
            config.pipeline.queue_memory_limit,
            1024 * 1024 * 1024,
            "1 GiB total → 1 GiB queue_memory_limit on the BamPipelineConfig",
        );
    }

    #[test]
    fn does_not_support_when_no_umi_source() {
        let dir = tempfile::TempDir::new().unwrap();
        // No UMIs supplied.
        let runall =
            runall_for_correct(dir.path().join("in.bam"), dir.path().join("out.bam"), Vec::new());
        let runner = CorrectChainRunner::from_runall(&runall);
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Correct)));
    }

    #[test]
    fn name_is_stable() {
        let dir = tempfile::TempDir::new().unwrap();
        let runall = runall_for_correct(
            dir.path().join("in.bam"),
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        let runner = CorrectChainRunner::from_runall(&runall);
        assert_eq!(runner.name(), CORRECT_CHAIN_RUNNER_NAME);
    }

    #[test]
    fn load_umi_sequences_from_strings_dedups_and_uppercases() {
        let umis = vec!["AAAA".to_string(), "CCCC".to_string(), "aaaa".to_string()];
        let (seqs, len) = load_umi_sequences(&umis, &[]).expect("should load UMIs");
        assert_eq!(len, 4);
        assert_eq!(seqs.len(), 2);
    }

    #[test]
    fn load_umi_sequences_empty_fails() {
        let result = load_umi_sequences(&[], &[]);
        assert!(result.is_err());
    }
}

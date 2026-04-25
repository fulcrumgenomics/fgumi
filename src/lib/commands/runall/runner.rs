//! Plan runner: instantiates a [`Plan`] into concrete pipeline types and runs it.
//!
//! [`run_plan`] dispatches on the combination of [`SourceSpec`] and [`SinkSpec`]
//! variants (2 × 3 = 6 pairs), builds the concrete [`Source`] and [`Sink`] for
//! each, converts every [`StageSpec`] into a [`StageEntry`], and invokes
//! [`run_multi_stage`] with the correctly typed `SourceOut`/`SinkIn` generics.
//!
//! This is the single owner of the "how" of running a pipeline. All callers
//! (runall entry points, tests, tools) just build a [`Plan`] and call
//! [`Plan::run`].
//!
//! [`Source`]: crate::runall::engine::source::Source
//! [`Sink`]: crate::runall::engine::sink::Sink

use std::path::PathBuf;

use anyhow::{Context, Result};
use fgumi_consensus::caller::ConsensusOutput;

use super::plan::{Plan, SinkSpec, SourceSpec, StageSpec};
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::deadlock::WatchdogConfig;
use crate::runall::engine::driver::{
    ErasedStage, StageEntry, TypedStage, run_multi_stage_with_metrics_and_watchdog,
};
use crate::runall::engine::fastq_types::PerStreamChunk;
use crate::runall::engine::output_types::{CompressedBatch, SerializedBatch, SerializedFastqBatch};
use crate::runall::engine::sink::Sink;
use crate::runall::engine::sink::bam::BamSink;
use crate::runall::engine::sink::bam_file_write::BamFileWrite;
use crate::runall::engine::sink::fastq_file::FastqFileSink;
use crate::runall::engine::source::fastq::FastqFileRead;
use crate::runall::engine::source::{BamBatchedSource, Source};
use crate::runall::engine::special_stage::TypedSpecialStage;
use crate::runall::engine::stage_source::{ErasedStageSource, StageSource};
use crate::runall::engine::stages::align_and_merge::AlignAndMerge;
use crate::runall::engine::stages::coalesce::Coalesce;
use crate::runall::engine::stages::group_assign::GroupAssignStage;
use crate::runall::engine::stages::position_batch::PositionBatchStage;
use crate::runall::engine::stages::sort::SortStage;
use crate::runall::engine::stages::{
    BgzfCompress, BgzfDecompress, ConsensusStage, CorrectStage, ExtractStage, FastqBlockMerge,
    FastqBlockParse, FastqPair, FastqParse, FilterStage, MiGroupBatch, ToFastq,
};

impl Plan {
    /// Instantiate and run this plan to completion.
    ///
    /// Blocks until the pipeline finishes (success, cancel, or error). The
    /// caller is responsible for any post-run finalization (metrics, temp →
    /// final rename, etc.).
    ///
    /// # Errors
    ///
    /// Returns an error if any stage fails or if a worker thread panics.
    #[tracing::instrument(name = "plan_run", skip_all, fields(num_stages = self.stages.len()))]
    pub fn run(self, cancel: &CancelToken) -> Result<()> {
        run_plan(self, cancel)
    }
}

/// Free-function equivalent of [`Plan::run`].
///
/// # Errors
///
/// Returns an error if any stage fails or if a worker thread panics.
pub fn run_plan(plan: Plan, cancel: &CancelToken) -> Result<()> {
    tracing::info!("pipeline plan: {}", plan.summary());

    let Plan { source, stages, sink, config, metrics_tsv } = plan;
    let entries: Vec<StageEntry> = stages.into_iter().map(spec_to_entry).collect();

    match (source, sink) {
        (SourceSpec::FastqFileRead { paths }, sink) => {
            let src: Box<dyn Source<Output = PerStreamChunk>> =
                Box::new(FastqFileRead::new(paths)?);
            run_with_source::<PerStreamChunk>(src, entries, sink, cancel, config, metrics_tsv)
        }
        (SourceSpec::BamBatchedSource { path, threads }, sink) => {
            let src: Box<dyn Source<Output = SerializedBatch>> =
                Box::new(BamBatchedSource::new(path, threads));
            run_with_source::<SerializedBatch>(src, entries, sink, cancel, config, metrics_tsv)
        }
    }
}

/// Generous no-progress timeout for the deadlock watchdog. Must cover the
/// slowest legitimate "no progress" interval in any supported pipeline —
/// today that's the initial bwa reference load, which can be >60s on cold
/// I/O against a WGS index. The watchdog only trips if *no* stage
/// makes any progress for this long, so a large value is safe.
fn runall_watchdog_cfg() -> WatchdogConfig {
    WatchdogConfig { timeout: std::time::Duration::from_secs(120), ..WatchdogConfig::default() }
}

/// Dispatch on [`SinkSpec`] once the source side has been concretized.
fn run_with_source<SourceOut>(
    source: Box<dyn Source<Output = SourceOut>>,
    entries: Vec<StageEntry>,
    sink: SinkSpec,
    cancel: &CancelToken,
    config: crate::runall::engine::driver::PipelineConfig,
    metrics_tsv: Option<PathBuf>,
) -> Result<()>
where
    SourceOut: Send + 'static,
{
    let wd = runall_watchdog_cfg();
    match sink {
        SinkSpec::BamFileWrite {
            path,
            header,
            secondary_path,
            secondary_header,
            queue_capacity,
        } => {
            let sink: Box<dyn Sink<Input = CompressedBatch>> = Box::new(BamFileWrite::new(
                path,
                header,
                secondary_path,
                secondary_header,
                queue_capacity,
            ));
            run_multi_stage_with_metrics_and_watchdog::<SourceOut, CompressedBatch>(
                source,
                entries,
                sink,
                cancel,
                config,
                metrics_tsv.as_deref(),
                wd,
            )
            .context("pipeline plan failed")
        }
        SinkSpec::FastqFileSink { path, queue_capacity } => {
            let sink: Box<dyn Sink<Input = SerializedFastqBatch>> =
                Box::new(FastqFileSink::new(path, queue_capacity));
            run_multi_stage_with_metrics_and_watchdog::<SourceOut, SerializedFastqBatch>(
                source,
                entries,
                sink,
                cancel,
                config,
                metrics_tsv.as_deref(),
                wd,
            )
            .context("pipeline plan failed")
        }
        SinkSpec::BamSink { path, header, threads, compression_level, deferred_header } => {
            let mut bam_sink = BamSink::new(path, header, threads, compression_level);
            if let Some(slot) = deferred_header {
                bam_sink = bam_sink.with_deferred_header(slot);
            }
            let sink: Box<dyn Sink<Input = ConsensusOutput>> = Box::new(bam_sink);
            run_multi_stage_with_metrics_and_watchdog::<SourceOut, ConsensusOutput>(
                source,
                entries,
                sink,
                cancel,
                config,
                metrics_tsv.as_deref(),
                wd,
            )
            .context("pipeline plan failed")
        }
    }
}

/// Convert a single [`StageSpec`] into a [`StageEntry`] the engine can dispatch.
fn spec_to_entry(spec: StageSpec) -> StageEntry {
    match spec {
        StageSpec::BgzfDecompress => pool_clone(BgzfDecompress::new()),
        StageSpec::BgzfCompress { level } => pool_factory(move || BgzfCompress::new(level)),
        StageSpec::FastqParse => pool_clone(FastqParse),
        StageSpec::FastqPair { num_streams } => pool_clone(FastqPair::new(num_streams)),
        StageSpec::FastqBlockParse => pool_clone(FastqBlockParse),
        StageSpec::FastqBlockMerge { num_streams } => pool_clone(FastqBlockMerge::new(num_streams)),
        StageSpec::Extract { params, read_structures, encoding } => {
            pool_clone(ExtractStage::new(params, read_structures, encoding))
        }
        StageSpec::Correct { state } => {
            let encoded = state.encoded.clone();
            let metrics = state.metrics.clone();
            let umi_length = state.umi_length;
            let max_mismatches = state.max_mismatches;
            let min_distance_diff = state.min_distance_diff;
            let revcomp = state.revcomp;
            let dont_store = state.dont_store_original_umis;
            let track_rejects = state.track_rejects;
            let cache_size = state.cache_size;
            pool_factory(move || {
                CorrectStage::new(
                    encoded.clone(),
                    umi_length,
                    max_mismatches,
                    min_distance_diff,
                    revcomp,
                    dont_store,
                    track_rejects,
                    cache_size,
                    metrics.clone(),
                )
            })
        }
        StageSpec::ToFastq => pool_clone(ToFastq::with_defaults(false)),
        StageSpec::AlignAndMerge {
            aligner_command,
            reference,
            dict_path,
            unmapped_header,
            tag_info,
            skip_pa_tags,
            no_read_suffix,
            deferred_header,
        } => {
            let mut stage = AlignAndMerge::new(
                aligner_command,
                reference,
                dict_path,
                unmapped_header,
                tag_info,
                skip_pa_tags,
                no_read_suffix,
            );
            if let Some(slot) = deferred_header {
                stage = stage.with_deferred_header(slot);
            }
            StageEntry::Special(Box::new(TypedSpecialStage::new(stage)))
        }
        StageSpec::Sort { memory_limit, temp_dir, threads, header } => {
            StageEntry::Special(Box::new(TypedSpecialStage::new(SortStage::new(
                memory_limit,
                temp_dir,
                threads,
                header,
            ))))
        }
        StageSpec::PositionBatch { header } => {
            StageEntry::Special(Box::new(TypedSpecialStage::new(PositionBatchStage::new(header))))
        }
        StageSpec::GroupAssign { umi_tag, assign_tag, strategy, max_edits, filter_config } => {
            pool_clone(GroupAssignStage::new(
                strategy,
                max_edits,
                umi_tag,
                assign_tag,
                filter_config,
            ))
        }
        StageSpec::Consensus { factory, call_overlapping_bases } => {
            pool_factory(move || ConsensusStage::with_overlapping(&factory, call_overlapping_bases))
        }
        StageSpec::Filter { config, filter_by_template, metrics } => {
            let mut stage = FilterStage::new(config, filter_by_template);
            if let Some(acc) = metrics {
                stage = stage.with_metrics(acc);
            }
            pool_clone(stage)
        }
        StageSpec::CoalesceSerialized { target_chunk_size } => {
            let c: Coalesce<SerializedBatch> =
                Coalesce::builder().target_chunk_size(target_chunk_size).build();
            StageEntry::Special(Box::new(TypedSpecialStage::new(c)))
        }
        StageSpec::CoalesceMiGroup { target_chunk_size } => {
            let c: Coalesce<MiGroupBatch> =
                Coalesce::builder().target_chunk_size(target_chunk_size).build();
            StageEntry::Special(Box::new(TypedSpecialStage::new(c)))
        }
    }
}

/// Wrap a `Clone` stage as a pool entry.
fn pool_clone<S>(stage: S) -> StageEntry
where
    S: crate::runall::engine::stage::Stage + Clone + Send + Sync + 'static,
{
    StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(stage)))
}

/// Wrap a factory as a pool entry (for stages that can't be `Clone`).
///
/// Mirrors `PipelineBuilder::stage_factory`: probes the factory once to
/// capture declared parallelism and stage name, then hands the factory to
/// [`ErasedStageSource::from_erased_factory`] so fresh instances are built
/// per worker.
fn pool_factory<S, F>(factory: F) -> StageEntry
where
    S: crate::runall::engine::stage::Stage + Send + 'static,
    F: Fn() -> S + Send + Sync + 'static,
{
    let probe = factory();
    let parallelism = probe.parallelism();
    let name = probe.name();
    drop(probe);
    StageEntry::Pool(ErasedStageSource::from_erased_factory(
        move || Box::new(TypedStage::new(factory())) as Box<dyn ErasedStage>,
        parallelism,
        name,
    ))
}

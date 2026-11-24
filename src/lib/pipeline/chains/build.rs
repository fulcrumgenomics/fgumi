//! [`build_for`] — the single chain-construction entry point.
//!
//! `build_for` constructs the chain as a stage-by-stage loop:
//!
//!   1. Validates the spec (progression, options presence, cross-stage constraints).
//!   2. Constructs a [`ChainBuilder`] from the spec.
//!   3. Calls `add_source`.
//!   4. Walks `spec.stages`, calling `add_stage(stage, position)` for each,
//!      where `position` is `Terminal` for the last stage and `Intermediate` for all
//!      earlier ones.
//!   5. Calls `add_sink`.
//!   6. Calls `build()` to produce the [`BuiltPipeline`]; its last step creates
//!      the rejects files the stages queued.
//!
//! Output files are created only when the build gets that far: rejects files
//! are deferred to the end of `build()`, and if any step fails, every output
//! file the build already created (the main output, any rejects) is removed
//! before the error is returned, so a failed build leaves nothing behind.
//!
//! Every chain — including the sole-`[Stage::Sort]` standalone sort — runs
//! through this uniform `add_source` → stages → `add_sink` topology on the
//! work-stealing pool. There is no longer a Sort-terminal special case: the
//! former self-contained file→file `SortBamFile` step was removed, and
//! standalone sort now streams (`ParseBamRecords` → arena sort ingest →
//! `SpillGather` → `SpillBlockCompress` → `SpillWrite` → `SortSpillDecompress` →
//! `SortMerge` → sink) like every other command.

use anyhow::{Result, bail};

use crate::pipeline::chains::{
    BuiltPipeline, ChainSpec,
    builder::{ChainBuilder, StagePosition},
    validate::{
        validate_cross_stage_constraints, validate_stage_opts_present, validate_stage_progression,
    },
};

/// Construct a [`BuiltPipeline`] from a [`ChainSpec`]. See the design
/// doc (`docs/design/refactor/2026-05-20-unified-chain-builder-design.md`)
/// for the architectural rationale.
///
/// Does NOT execute the pipeline — the caller picks `PipelineConfig`
/// (threads, scheduler) and calls `built.pipeline.run(built.config)`,
/// then drains `built.finalize` in order. The convenience method
/// [`BuiltPipeline::run`] wraps both steps.
///
/// ## Dispatch
///
/// `build_for` constructs a [`ChainBuilder`] and drives it stage-by-stage:
///
/// ```text
/// chain.add_source()?
/// for (i, stage) in stages.iter().enumerate() {
///     let pos = if i == last { Terminal } else { Intermediate };
///     chain.add_stage(*stage, pos)?;
/// }
/// chain.add_sink()?
/// chain.build()
/// ```
///
/// No output file survives a build that returns `Err` (see the module docs).
///
/// Each `add_<stage>` method reads its options from `spec.stage_opts` and pushes
/// the appropriate typed-step sequence onto the pipeline builder. The position
/// argument controls whether a serialise-to-bytes step is appended (Terminal)
/// or whether the native typed output is left for the next stage (Intermediate).
///
/// # Errors
///
/// Returns the first validator violation, or any underlying step
/// construction error (header parse, file-open, etc.).
pub fn build_for(spec: ChainSpec) -> Result<BuiltPipeline> {
    build_for_with_source_validation(spec, |_| Ok(()))
}

/// Build a chain after validating its detected source format. The source has
/// been opened once and retains its prefix; no stage or sink output is opened
/// before `validate_source` succeeds.
///
/// # Errors
///
/// Source validation or any normal [`build_for`] construction failure.
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn build_for_with_source_validation<F>(
    spec: ChainSpec,
    validate_source: F,
) -> Result<BuiltPipeline>
where
    F: FnOnce(bool) -> Result<()>,
{
    validate_stage_progression(&spec)?;
    validate_stage_opts_present(&spec)?;
    validate_cross_stage_constraints(&spec)?;
    // Re-check source invariants the variant types cannot encode (e.g. a
    // struct-literal `Fastqs` with mismatched paths/read_structures lengths),
    // so the safety net documented on `SourceSpec::validate` is delivered on
    // the build path rather than left to construction-site discipline.
    spec.source.validate()?;

    // Defensive: `validate_stage_progression` above already rejects an empty
    // `stages`, so this guard is unreachable in practice. It is kept to
    // protect the `spec.stages.len() - 1` below from underflowing if the
    // validator order is ever changed.
    if spec.stages.is_empty() {
        bail!("ChainSpec must specify at least one stage");
    }

    let mut chain = ChainBuilder::new(&spec)?;
    validate_source(chain.source_is_cram())?;

    // NOTE: `--threads 1` whole-chain fusion is NOT wired here. It is applied
    // generically by the runtime (`Pipeline::run` → `run_fused_single_thread`)
    // to ANY linear source→sink chain at one worker — see
    // `run_fused_single_thread` in fgumi-pipeline-core. `build_for` therefore constructs the
    // same staged chain regardless of thread count; the runtime fuses it.
    //
    // Every chain — including the sole-`[Stage::Sort]` chain — runs through the
    // normal `add_source` → stages → `add_sink` topology on the work-stealing
    // pool. The former Sort-terminal special case (which skipped source/sink so
    // a self-contained `SortBamFile` could read/write files itself) was removed:
    // standalone sort now streams like every other command.
    chain.add_source()?;

    let last_idx = spec.stages.len() - 1;
    for (i, &stage) in spec.stages.iter().enumerate() {
        let position =
            if i == last_idx { StagePosition::Terminal } else { StagePosition::Intermediate };
        chain.add_stage(stage, position)?;
    }

    chain.add_sink()?;

    chain.build()
}

/// What [`build_merge`] builds: a k-way merge of sorted BAM inputs
/// (`fgumi merge`).
pub struct MergeSpec {
    /// The inputs, in input order (the merge's tie-break order), each opened
    /// once by the caller (a stream's header already parsed through its
    /// source).
    pub inputs: Vec<fgumi_pipeline_io::sort::MergeInput>,
    /// The output header, written as is.
    pub header: noodles::sam::Header,
    /// The order the inputs are sorted in.
    pub order: fgumi_sort::SortOrder,
    /// That order's `--order` value, for the unsorted-input error's hint
    /// ([`SortOrderArg::flag_value`](crate::commands::sort::SortOrderArg::flag_value)).
    pub order_name: &'static str,
    /// The full template key's context (template-coordinate only), from the
    /// merged header.
    pub template_ctx: Option<std::sync::Arc<fgumi_sort::TemplateKeyContext>>,
    /// A rewrite of each record before its key is extracted (read-tag
    /// renames), if any input needs one.
    pub rewrite: Option<fgumi_sort::RecordRewrite>,
    /// Where the merged alignment file is written (a staged temp, or stdout).
    pub output: std::path::PathBuf,
    /// The user's destination, used to select the format and name errors.
    pub output_name: std::path::PathBuf,
    /// Output BGZF compression level.
    pub compression: u32,
    /// `--threads`: how many threads compress the output. The pool adds up
    /// to as many workers again, one per input, so reading and decompression
    /// overlap output compression even at one thread (see
    /// `merge_pool_threads`).
    pub threads: usize,
    /// `--deadlock-timeout`: the deadlock watchdog's timeout, in seconds
    /// (`0` disables it).
    pub deadlock_timeout_secs: u64,
    /// Filled with the merge's record counts on success.
    pub stats_slot: std::sync::Arc<parking_lot::Mutex<Option<fgumi_sort::SortStats>>>,
}

/// Read-ahead budget per `--threads`: what `fgumi sort` spends per thread by
/// default, so a merge reads ahead as the sort's own merge does. Scaled by
/// `--threads`, not by the pool's worker count: the extra workers read what
/// this budget lets the planner request, they do not need more of it.
const MERGE_MEMORY_PER_THREAD: usize = 768 << 20;

/// Inter-step queue budget per merge pool worker. Not the sort's 768 MiB: the
/// merge is bound by output compression, so a large queue budget (and the
/// rebalancer feeding the fullest edge) parks merged batches in front of
/// `BgzfCompress` without making it any faster. The total is shared by every
/// byte-bounded queue, each floored at 1 MiB, so the queues can hold slightly
/// more than it at one or two workers.
///
/// 16 MiB is a choice, not a tuned optimum. On one local three-input merge
/// (three ~1.8 GB coordinate inputs, `--compression-level 1`, warm cache, one
/// run per value, a build of this chain with a single extra pool worker),
/// per-worker budgets of 4, 8, 16, 32 and 64 MiB gave wall times within 0.5%
/// of each other at `--threads` 1 and 8, while peak RSS at `--threads 8` grew
/// with the budget (162 MB at 4 MiB to 280 MB at 64 MiB). Splitting one of
/// those inputs into 50 (`--threads` + min(50, `--threads`) workers, two runs
/// each), 4, 16 and 64 MiB gave wall times within 5% of each other at
/// `--threads` 1 and 8. Independent of [`MERGE_MEMORY_PER_THREAD`], which
/// still sizes the read-ahead.
const MERGE_QUEUE_MEMORY_PER_THREAD: usize = 16 << 20;

/// The merge's pool size for `threads` (`--threads`) compression threads over
/// `inputs` inputs: one extra worker per input, at most `threads` and at least
/// one. `BgzfCompress` is capped at `threads`, so at most `threads` workers
/// compress at once and the rest read and decompress. The engine this chain
/// replaced read and decompressed every input on its own thread beside the
/// `threads` compressors; without an extra worker, `--threads 1` ran read,
/// decompress and compress back to back on one worker, and `--threads N` had
/// N workers do the work of N + inputs threads.
#[must_use]
fn merge_pool_threads(threads: usize, inputs: usize) -> usize {
    let threads = threads.max(1);
    threads + inputs.clamp(1, threads)
}

/// Attach the merge's format-specific sink, encoding on the existing worker pool.
fn add_merge_sink(
    tail: fgumi_pipeline_core::builder::Chain<
        '_,
        fgumi_pipeline_core::Single<fgumi_pipeline_io::types::DecompressedBlock>,
    >,
    output: &std::path::Path,
    output_name: &std::path::Path,
    header: &noodles::sam::Header,
    compression: u32,
    threads: usize,
    byte_limit: u64,
) -> Result<()> {
    use crate::pipeline::steps::bgzf::compress::BgzfCompress;
    use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;
    #[cfg(feature = "cram")]
    use std::sync::Arc;

    // At most `threads` workers encode; the pool's extra ones read and decompress.
    let cap = Some(fgumi_pipeline_core::PhaseCap::new("merge-compress", threads));
    if fgumi_bam_io::is_cram_path(output_name) {
        #[cfg(feature = "cram")]
        {
            use fgumi_pipeline_io::sink::write_cram::{
                CramEncode, CramHeader, CramPack, CramSetup, CramWrite, FinisherSlot,
            };
            use fgumi_pipeline_io::types::DecompressedBlock;

            let policy = fgumi_bam_io::CramOutputPolicy::new(
                output_name,
                !header.reference_sequences().is_empty(),
                compression,
            );
            let setup = CramSetup {
                path: output.to_path_buf(),
                display: output_name.to_path_buf(),
                header: CramHeader::Eager(Box::new(header.clone())),
                reference: policy.reference,
                no_ref: policy.no_ref,
                file_id: policy.file_id,
                crai: None,
            };
            let finisher: FinisherSlot = Arc::default();
            tail.chain(CramPack::<DecompressedBlock>::new(
                setup,
                Arc::clone(&finisher),
                byte_limit,
            ))
            .chain(CramEncode::new(byte_limit).with_phase_cap(cap))
            .chain(CramWrite::new(finisher))
            .into_sink_marker();
        }
    } else {
        let compress = BgzfCompress::new(compression, byte_limit, false).with_phase_cap(cap);
        let write = WriteBgzfFile::new(output, header, compression)
            .map_err(|e| anyhow::anyhow!("WriteBgzfFile::new: {e}"))?;
        tail.chain(compress).chain(write).into_sink_marker();
    }
    Ok(())
}

/// Build `fgumi merge`'s chain: the sort's phase 2 over the inputs.
///
/// ```text
/// AnnounceMergeInputs ─▶ SpillReadPlanner ─┬─ b0 ─▶ PreadSpillSlices ─▶ SortSpillDecompress
///                                          ├─ b1 ─▶ SortMerge<BlockOutput> ─▶ BgzfCompress ─▶ WriteBgzfFile
///                                          └─ b2 ─▶ NoConsolidation
/// ```
///
/// Each input is a body-keyed slot (its header skipped, its records keyed from
/// their bodies); the merge checks every input is sorted in `order` and merges
/// template-coordinate inputs on the full key. There is no consolidation: the
/// merge opens every input at once.
///
/// # Errors
///
/// The output cannot be created, or the pipeline does not build.
pub fn build_merge(spec: MergeSpec) -> Result<BuiltPipeline> {
    use std::sync::Arc;

    use crate::pipeline::steps::tuning::BamPipelineTuning;
    use fgumi_pipeline_io::sort::{
        AnnounceMergeInputs, BlockOutput, MergeInputs, NoConsolidation, SortMerge,
        SortSpillDecompress, SpillReadPlanner, SpillSupply, serial_merge_output,
    };

    let MergeSpec {
        inputs,
        header,
        order,
        order_name,
        template_ctx,
        rewrite,
        output,
        output_name,
        compression,
        threads,
        deadlock_timeout_secs,
        stats_slot,
    } = spec;
    fgumi_bam_io::validate_alignment_output(&output_name)?;
    let threads = threads.max(1);
    let pool_threads = merge_pool_threads(threads, inputs.len());
    let byte_limit = BamPipelineTuning::auto_tuned(threads).per_step_byte_limit;
    let names: Vec<String> = inputs.iter().map(|i| i.path.display().to_string()).collect();
    let k = inputs.len();
    let supply = SpillSupply::new();
    // Every pool worker can read and decompress (only compression is capped),
    // so the slice pool and the planner's in-flight cap (twice its thread
    // count) are sized from the pool, as the sort sizes them from the workers
    // that run its phase 2. Sized from `threads`, `--threads 1` would allow two
    // slices in flight across two reading workers: no slice queued behind the
    // one each worker is reading.
    let slices = fgumi_bam_io::pread::SliceBufferPool::new(2 * pool_threads + 8);
    let eligible = fgumi_pipeline_core::runtime::parallel_hosts(
        fgumi_pipeline_core::PoolPlacement::AllWorkers,
        None,
        pool_threads,
        None,
    )
    .clone_count();
    let total_memory = (MERGE_MEMORY_PER_THREAD * threads) as u64;
    let planner = SpillReadPlanner::new(
        total_memory,
        pool_threads,
        eligible,
        byte_limit,
        &supply,
        Arc::clone(&slices),
    )
    .with_read_streams(fgumi_bam_io::pread::ReadStreamsPolicy::auto());
    let pread = fgumi_pipeline_io::pread::PreadSlices::spill(slices, byte_limit);
    let decompress = SortSpillDecompress::new(byte_limit, &supply);
    let merge = SortMerge::<BlockOutput>::new(order, byte_limit)
        .with_spill_supply(&supply)
        .with_stats_slot(stats_slot)
        .with_merge_inputs(MergeInputs { inputs: names, order_name, template_ctx, rewrite });
    let builder = fgumi_pipeline_core::builder::PipelineBuilder::new();
    let fan =
        builder.chain(AnnounceMergeInputs::new(inputs, byte_limit)).chain(planner).into_multi();
    fan.b0.chain(pread).chain(decompress).into_sink_marker();
    // `fgumi merge` never takes the single-chunk fast path (its inputs are
    // slots), so the merge's job branch gets an idle `MergeWork`.
    let tail = serial_merge_output(fan.b1.chain(merge));
    add_merge_sink(tail, &output, &output_name, &header, compression, threads, byte_limit)?;
    fan.b2.chain(NoConsolidation).into_sink_marker();
    let pipeline = builder.build().map_err(|e| anyhow::anyhow!("Pipeline::build: {e}"))?;
    // The runtime every sort-chain command shares (the deadlock monitor on),
    // with the merge's own small queue budget rather than the sort's
    // `--max-memory` default (merge has no memory flag).
    let (scheduler, queue_memory) = crate::commands::common::sort_chain_runtime(
        crate::commands::common::MemoryLimit::Fixed(MERGE_QUEUE_MEMORY_PER_THREAD),
        crate::commands::common::MemoryReserve::Auto,
        true,
    );
    let scheduler = crate::commands::common::SchedulerOptions {
        deadlock: crate::commands::common::DeadlockTimeoutOptions {
            deadlock_timeout: deadlock_timeout_secs,
        },
        ..scheduler
    };
    let (config, _stats) = crate::pipeline::chains::build_helpers::build_pipeline_config_for_chain(
        &pipeline,
        pool_threads,
        &scheduler,
        &queue_memory,
    )?;
    log::info!(
        "Threads: {threads} (output compression) + {} reading and decompressing the {k} \
         inputs = {pool_threads} pool workers",
        pool_threads - threads
    );
    Ok(BuiltPipeline {
        pipeline,
        config,
        finalize: Vec::new(),
        finalize_on_success: Vec::new(),
        sort_stats_slot: None,
        #[cfg(test)]
        merge_class_tables: Vec::new(),
        #[cfg(test)]
        consolidation_class_table: None,
        #[cfg(test)]
        sort_caps: [None, None],
        #[cfg(test)]
        spill_eligible_clones: Some(eligible),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A merge spec over `k` positional inputs. The chain only stores its
    /// inputs while it is built, so each is an arbitrary non-empty file.
    fn merge_spec(dir: &std::path::Path, threads: usize, k: usize) -> MergeSpec {
        let inputs = (0..k)
            .map(|i| {
                let path = dir.join(format!("in{i}.bam"));
                std::fs::write(&path, [0u8; 64]).expect("write input");
                let file = std::fs::File::open(&path).expect("open input");
                fgumi_pipeline_io::sort::MergeInput {
                    path,
                    source: std::sync::Arc::new(file),
                    len: Some(64),
                    prefix_len: None,
                }
            })
            .collect();
        MergeSpec {
            inputs,
            header: noodles::sam::Header::default(),
            order: fgumi_sort::SortOrder::Coordinate,
            order_name: crate::commands::sort::SortOrderArg::Coordinate.flag_value(),
            template_ctx: None,
            rewrite: None,
            output: dir.join("out.bam"),
            output_name: dir.join("out.bam"),
            compression: 1,
            threads,
            deadlock_timeout_secs: crate::commands::common::DEFAULT_DEADLOCK_TIMEOUT_SECS,
            stats_slot: std::sync::Arc::new(parking_lot::Mutex::new(None)),
        }
    }

    #[cfg(not(feature = "cram"))]
    #[test]
    fn merge_cram_refusal_names_user_destination_without_creating_output() {
        let dir = tempfile::tempdir().unwrap();
        let mut spec = merge_spec(dir.path(), 1, 2);
        spec.output = dir.path().join("staged.tmp");
        spec.output_name = dir.path().join("user.cram");
        let output = spec.output.clone();
        let named = spec.output_name.clone();
        let err = build_merge(spec).err().expect("CRAM requires the opt-in feature");
        assert!(err.to_string().contains(&named.display().to_string()));
        assert!(err.to_string().contains("--features cram"));
        assert!(!output.exists());
    }

    /// The merge emits CRAM directly, including the eager header, records,
    /// user-facing file ID and encoder cap, even when its write path is staged.
    #[cfg(feature = "cram")]
    #[rstest::rstest]
    fn merge_writes_cram(#[values(1, 4)] threads: usize, #[values(0, 12_050)] count: usize) {
        use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
        use noodles::sam::alignment::io::Write as _;
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        let mut spec = merge_spec(dir.path(), threads, 0);
        spec.header = noodles::sam::Header::builder().add_comment("merge header").build();
        let mut expected = Vec::new();
        for shard in 0..2 {
            let path = dir.path().join(format!("in{shard}.bam"));
            let mut writer = noodles::bam::io::Writer::new(std::fs::File::create(&path).unwrap());
            writer.write_header(&spec.header).unwrap();
            for i in (shard..count).step_by(2) {
                let record = SamBuilder::new()
                    .read_name(format!("read{i:07}").as_bytes())
                    .flags(flags::UNMAPPED)
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8])
                    .build();
                let buf = fgumi_raw_bam::raw_record_to_record_buf(&record, &spec.header).unwrap();
                writer.write_alignment_record(&spec.header, &buf).unwrap();
                expected.push(record);
            }
            writer.try_finish().unwrap();
            drop(writer);
            let file = std::fs::File::open(&path).unwrap();
            let len = file.metadata().unwrap().len();
            spec.inputs.push(fgumi_pipeline_io::sort::MergeInput {
                path,
                source: Arc::new(file),
                len: Some(len),
                prefix_len: None,
            });
        }
        spec.order = fgumi_sort::SortOrder::Queryname(fgumi_sort::QuerynameComparator::Natural);
        spec.order_name = "queryname";
        spec.output = dir.path().join("staged.tmp");
        spec.output_name = dir.path().join("merged.cram");
        let output = spec.output.clone();
        let header = spec.header.clone();
        let built = build_merge(spec).unwrap();
        let caps = built.pipeline.phase_caps();
        assert_eq!(caps.len(), 1);
        assert_eq!(caps[0].0, "CramEncode");
        assert_eq!(caps[0].1.max(), threads);
        built.run().expect("CRAM merge");
        let bytes = std::fs::read(&output).unwrap();
        assert_eq!(&bytes[..4], b"CRAM");
        assert!(bytes[6..26].starts_with(b"merged.cram"));
        let (mut reader, actual_header) = fgumi_bam_io::create_raw_bam_reader_with_opts(
            &output,
            1,
            fgumi_bam_io::PipelineReaderOpts::default(),
        )
        .unwrap();
        assert_eq!(actual_header.comments(), header.comments());
        expected.sort_by(|a, b| a.read_name().cmp(b.read_name()));
        let mut record = RawRecord::default();
        for expected in expected {
            assert_ne!(reader.read_record(&mut record).unwrap(), 0);
            assert_eq!(record.read_name(), expected.read_name());
            assert_eq!(record.flags(), expected.flags());
            assert_eq!(record.sequence_vec(), expected.sequence_vec());
            assert_eq!(record.quality_scores(), expected.quality_scores());
        }
        assert_eq!(reader.read_record(&mut record).unwrap(), 0);
    }

    /// `fgumi merge --deadlock-timeout` reaches the merge pipeline's watchdog
    /// (issue #1068: it used to be hardcoded).
    #[rstest::rstest]
    #[case::default(crate::commands::common::DEFAULT_DEADLOCK_TIMEOUT_SECS)]
    #[case::raised(120)]
    #[case::disabled(0)]
    fn build_merge_honours_the_deadlock_timeout(#[case] secs: u64) {
        let dir = tempfile::tempdir().expect("temp dir");
        let spec = MergeSpec { deadlock_timeout_secs: secs, ..merge_spec(dir.path(), 2, 2) };
        let built = build_merge(spec).expect("build_merge");
        assert_eq!(built.config.deadlock_timeout_secs, secs);
    }

    /// The pool runs `--threads` plus one worker per input (at most
    /// `--threads` more); only `BgzfCompress` is capped, at `--threads`, so
    /// the extra workers read and decompress; and the queue budget is the
    /// merge's own small per-worker budget, not the read-ahead's.
    #[rstest::rstest]
    #[case::one_thread_two_inputs(1, 2, 2)]
    #[case::eight_threads_two_inputs(8, 2, 10)]
    #[case::eight_threads_twelve_inputs(8, 12, 16)]
    fn build_merge_sizes_its_pool_cap_and_queue_budget(
        #[case] threads: usize,
        #[case] k: usize,
        #[case] pool: usize,
    ) {
        let dir = tempfile::tempdir().expect("temp dir");
        let built = build_merge(merge_spec(dir.path(), threads, k)).expect("build_merge");
        assert_eq!(merge_pool_threads(threads, k), pool);
        assert_eq!(built.config.threads, pool, "--threads plus one worker per input");
        let caps = built.pipeline.phase_caps();
        let capped: Vec<&str> = caps.iter().map(|(name, _)| *name).collect();
        assert_eq!(capped, ["BgzfCompress"], "only output compression is capped");
        assert_eq!(caps[0].1.max(), threads);
        let budget = built.config.queue_memory_total.expect("a queue budget");
        assert_eq!(budget, (MERGE_QUEUE_MEMORY_PER_THREAD * pool) as u64);
        assert!(budget < (MERGE_MEMORY_PER_THREAD * threads) as u64, "not the read-ahead's budget");
    }

    /// The planner's slice-per-fill bound (`eligible`) must equal the built
    /// pipeline's real placement plan for `PreadSpillSlices`: the pool's
    /// size, and so the clone count, is decided separately by
    /// `build_pipeline_config_for_chain`. As
    /// `spill_eligible_clones_match_the_planned_placement` pins the sort's.
    #[rstest::rstest]
    #[case::one_thread(1, 2)]
    #[case::eight_threads(8, 2)]
    #[case::eight_threads_many_inputs(8, 12)]
    fn spill_eligible_clones_match_the_planned_placement(#[case] threads: usize, #[case] k: usize) {
        let dir = tempfile::tempdir().expect("temp dir");
        let built = build_merge(merge_spec(dir.path(), threads, k)).expect("build_merge");
        let planned = built.pipeline.planned_clone_counts(built.config.threads);
        let clones = planned
            .iter()
            .find(|(name, _)| *name == "PreadSpillSlices")
            .map(|&(_, c)| c)
            .expect("the input reads");
        assert_eq!(built.spill_eligible_clones, Some(clones), "{planned:?}");
    }
}

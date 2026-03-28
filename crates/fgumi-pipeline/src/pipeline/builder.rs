//! Pipeline builder: composes Zone 3 stages, wires queues, and runs the scheduler.
//!
//! Supports three entry points:
//! - `Stage::Group` — reads a template-coordinate-sorted BAM and feeds Zone 3 directly.
//! - `Stage::Sort` — reads an unsorted BAM, sorts by template-coordinate, then feeds Zone 3.
//! - `run_two_phase` — runs Phase A (extract → align → zipper → sort) then Phase B (Zone 3).
//!
//! All threads (including the calling thread during completion) participate in work-stealing.
//! There is no dedicated writer thread; writing is a stage runner with the highest priority.
//! This unified model works with any thread count, including `--threads 1`.

use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use anyhow::{Context, Result};
use log::info;

use fgumi_bgzf::reader::BGZF_EOF;
use fgumi_consensus::filter::FilterConfig;
use fgumi_lib::bam_io::create_raw_bam_reader;
use fgumi_lib::extract::{ExtractParams, QualityEncoding};
use fgumi_lib::grouper::GroupFilterConfig;
use fgumi_lib::sort::{
    LibraryLookup, RawExternalSorter, SortOrder, SortedRecordSink, extract_template_key_inline,
    merge_sorted_chunks_to_sink,
};
use fgumi_lib::template::TemplateIterator;
use fgumi_lib::unified_pipeline::MemoryEstimate;
use fgumi_umi::{TagInfo, UmiAssigner};

use crate::pipeline::batch_state::AlignerBatchState;
use crate::pipeline::phase_a::PhaseAStep;
use crate::pipeline::phase_a_graph::{ExtractSource, MergeFn, PhaseAGraph};
use crate::pipeline::phase_a_steps;
use crate::pipeline::phase_a_worker::{start_phase_a_pool, wait_phase_a_pool};
use crate::stages::aligner::AlignerProcess;

use crate::pipeline::queue::WorkQueue;
use crate::pipeline::reorder::ReorderBuffer;
use crate::pipeline::scheduler::{
    BackpressureState, PipelineStep, Scheduler, StageGraph, StepResult, WorkerState,
    WorkerStateSlot, WriteState,
};
use crate::pipeline::stats::StageStats;
use crate::stage::SequencedItem;
use crate::stages::compress::{CompressStage, CompressedBatch};
use crate::stages::consensus::{CallerFactory, ConsensusStage};
use crate::stages::filter::FilterStage;
use crate::stages::group_assign::{GroupAndAssignStage, MiGroup, PositionGroupBatch};
use crate::stages::position_batch::PositionBatcher;

// ============================================================================
// Configuration types
// ============================================================================

/// Callback type that builds the final output header from the mapped SAM header.
///
/// The stdout-reader thread calls this after reading the SAM header from the aligner.
/// The closure typically merges the unmapped header, mapped header, and reference dict,
/// then adds a PG record.
pub type BuildOutputHeaderFn =
    Box<dyn FnOnce(&noodles::sam::Header) -> anyhow::Result<noodles::sam::Header> + Send>;

/// Pipeline entry point: which stage to start from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage {
    /// Start from sorting: read unsorted BAM → template-coordinate sort → Zone 3.
    Sort,
    /// Start from grouping: read already-sorted BAM → Zone 3.
    Group,
}

/// Configuration for the runall pipeline.
pub struct PipelineConfig {
    /// Input BAM file path.
    pub input: PathBuf,
    /// Output BAM file path.
    pub output: PathBuf,
    /// Reference FASTA path (for consensus calling).
    pub reference: PathBuf,
    /// Total threads available to the pipeline.
    pub threads: usize,
    /// Which stage to begin from.
    pub start_from: Stage,

    // -- Stage configs --
    /// UMI assigner for grouping.
    pub assigner: Arc<dyn UmiAssigner>,
    /// Two-byte BAM tag containing the UMI (e.g., `b"RX"`).
    pub umi_tag: [u8; 2],
    /// Two-byte BAM tag for the assigned molecule ID (e.g., `b"MI"`).
    pub assign_tag: [u8; 2],
    /// Template filter config for pre-UMI-assignment filtering (MAPQ, unmapped, etc.).
    pub group_filter_config: GroupFilterConfig,
    /// Factory that creates per-thread consensus callers.
    pub caller_factory: CallerFactory,
    /// Filter configuration for consensus-read filtering.
    pub filter_config: FilterConfig,
    /// Minimum base quality for masking.
    pub min_base_quality: u8,
    /// BGZF compression level for output (1–12).
    pub compression_level: u32,

    // -- Sort config (used when start_from == Sort) --
    /// Memory limit for external sort (bytes).
    pub sort_memory_limit: usize,
    /// Temporary directory for sort spill files.
    pub sort_temp_dir: Option<PathBuf>,
}

/// Configuration for the full two-phase pipeline: Phase A (extract → align → zipper → sort)
/// followed by Phase B (Zone 3: group → consensus → filter → compress → write).
///
/// Phase A uses a work-stealing pool for Extract, `ToFastq`, and Zipper, with dedicated
/// I/O threads for the aligner stdout reader and sort-accumulate. Phase B reuses the
/// existing [`PipelineBuilder::run_from_sorted_chunks`] infrastructure.
///
/// The caller is responsible for:
/// 1. Spawning the aligner process and taking its stdin.
/// 2. Reading the SAM header from aligner stdout (to resolve the output header).
/// 3. Passing the already-header-consumed SAM reader and the mapped header to this config.
pub struct TwoPhaseConfig {
    // -- Phase A inputs --
    /// FASTQ source iterator for Extract.
    pub extract_source: ExtractSource,
    /// Extract parameters (tags, read group, read structures).
    pub extract_params: ExtractParams,
    /// Detected quality encoding for FASTQ input.
    pub quality_encoding: QualityEncoding,
    /// Aligner subprocess handle (stdin already taken; stdout already taken).
    /// Caller must pass this so `run_two_phase` can `wait()` on it.
    pub aligner_process: AlignerProcess,
    /// Aligner stdin pipe (taken from the aligner process before constructing this config).
    pub aligner_stdin: std::process::ChildStdin,
    /// Aligner stdout pipe (taken from the aligner process). The stdout-reader thread will
    /// read the SAM header and then records from this pipe.
    pub aligner_stdout: std::process::ChildStdout,
    /// Callback that builds the final output header from the mapped header read by the
    /// stdout-reader thread. This lets the binary crate inject its header-building logic
    /// (e.g. `build_output_header` + `add_pg_record`) without the library crate depending
    /// on it.
    pub build_output_header_fn: BuildOutputHeaderFn,
    /// Tag manipulation info for zipper merge.
    pub tag_info: Arc<TagInfo>,
    /// Whether to skip PA (primary alignment) tag generation in zipper.
    pub skip_pa_tags: bool,
    /// Zipper merge function: transfers tags from unmapped to mapped template.
    pub merge_fn: MergeFn,
    /// Batch coordination state for bwa `-K` flow control (optional).
    pub batch_state: Option<Arc<AlignerBatchState>>,

    // -- Correction config (optional) --
    /// UMI correction configuration. When `Some`, the Correct pool step is active
    /// between Extract and `ToFastq`, using per-worker LRU caches.
    pub correction_config: Option<crate::pipeline::phase_a_graph::CorrectionConfig>,

    // -- Sort config --
    /// Memory limit for external sort (bytes).
    pub sort_memory_limit: usize,
    /// Temporary directory for sort spill files.
    pub sort_temp_dir: Option<PathBuf>,
    /// Number of threads for the external sorter.
    pub sort_threads: usize,

    // -- Phase B config --
    /// Pipeline configuration for Zone 3 stages.
    pub pipeline_config: PipelineConfig,

    // -- Thread control --
    /// Total threads for the Phase A worker pool.
    pub threads: usize,
    /// Shared cancellation flag.
    pub cancel: Arc<AtomicBool>,
}

// ============================================================================
// Queue capacity defaults
// ============================================================================

/// Default channel capacity for inter-stage queues.
///
/// Sized large enough (4096) to absorb producer bursts without triggering the
/// spin-loop in [`WorkQueue::push_blocking`].  The real backpressure mechanism is the
/// per-queue memory limit, not the slot count.  At 5 queues the pre-allocated
/// overhead is ~5 * 4096 * 8 bytes ≈ 160 KB — negligible.
const QUEUE_CAPACITY: usize = 4096;

/// Default memory limit per queue (256 MB).
const QUEUE_MEMORY_LIMIT: usize = 256 * 1024 * 1024;

// ============================================================================
// PipelineBuilder
// ============================================================================

/// Builds and executes the runall pipeline from a [`PipelineConfig`].
///
/// The pipeline wires Zone 3 stages (group-assign → consensus → filter → compress → write)
/// with bounded work queues and runs them under a work-stealing scheduler. All threads,
/// including the calling thread during completion, participate in the same work-stealing
/// loop. Writing is the highest-priority stage runner, protected by a mutex.
pub struct PipelineBuilder;

impl PipelineBuilder {
    /// Build and run the pipeline according to `config`.
    ///
    /// # Errors
    ///
    /// Returns an error if any stage fails, I/O fails, or the scheduler encounters an error.
    pub fn run(config: &PipelineConfig) -> Result<()> {
        match config.start_from {
            Stage::Group => Self::run_from_group(config),
            Stage::Sort => Self::run_from_sort(config),
        }
    }

    /// Run the pipeline starting from an already template-coordinate-sorted BAM.
    ///
    /// Reads records from the input BAM, batches them by position into
    /// [`PositionGroupBatch`] items, and feeds them through Zone 3 stages.
    fn run_from_group(config: &PipelineConfig) -> Result<()> {
        info!("Pipeline: starting from grouped/sorted BAM input");

        // 1. Open input BAM and read header.
        // Single I/O thread for BAM reading; parallelism is in the worker pool.
        let io_threads = 1;
        let (reader, header) =
            create_raw_bam_reader(&config.input, io_threads).context("opening input BAM")?;
        let lib_lookup = LibraryLookup::from_header(&header);

        // 2. Build Zone 3 queues and stages, start scheduler (including write stage).
        let zone3 = Zone3::start(config, &header)?;

        // 3. Feed input records to the position batcher.
        let mut batcher = PositionBatcher::new(zone3.q_position_batch.clone());
        let read_ahead = fgumi_lib::sort::read_ahead::RawReadAheadReader::new(reader);

        for record in read_ahead {
            let bam_bytes = record.as_ref();
            let key = extract_template_key_inline(bam_bytes, &lib_lookup, None);
            batcher.emit(&key, bam_bytes.to_vec())?;
        }
        batcher.finish()?;
        let total_batches = batcher.batches_pushed();
        info!("Pipeline: fed {total_batches} position batches to Zone 3");

        // 4. Signal that no more position batches will arrive, then wait for completion.
        zone3.q_position_batch.close();
        zone3.wait_for_completion(total_batches)?;

        info!("Pipeline: done");
        Ok(())
    }

    /// Run the pipeline starting from an unsorted BAM (sort → Zone 3).
    ///
    /// Uses a two-phase approach:
    /// 1. **Sort-accumulate** (blocking on calling thread): reads the input BAM, sorts chunks
    ///    to disk/memory.
    /// 2. **Merge** (spawned on a feeder thread): K-way merges the sorted chunks through a
    ///    [`PositionBatcher`] into Zone 3.
    ///
    /// While the feeder thread drives the merge, the calling thread immediately joins the
    /// Zone 3 worker pool to help with consensus and other stages. This eliminates the
    /// blocking merge-on-calling-thread pattern where the calling thread was stuck on
    /// merge I/O instead of contributing to downstream work.
    fn run_from_sort(config: &PipelineConfig) -> Result<()> {
        info!("Pipeline: starting from unsorted BAM (sort first)");

        // 1. Open input to read the header for queue/stage setup.
        let (_, header) =
            create_raw_bam_reader(&config.input, 1).context("opening input BAM for header")?;

        // 2. Sort-accumulate phase (blocking on calling thread).
        let sort_threads = 2.max(config.threads.saturating_sub(2));
        let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(config.sort_memory_limit)
            .threads(sort_threads);
        if let Some(ref temp_dir) = config.sort_temp_dir {
            sorter = sorter.temp_dir(temp_dir.clone());
        }

        let sorted_chunks = sorter
            .sort_accumulate_template(&config.input, &header)
            .context("sort-accumulate phase")?;
        info!(
            "Pipeline: sort-accumulate complete ({} records, {} disk chunks)",
            sorted_chunks.stats.total_records, sorted_chunks.stats.chunks_written
        );

        // 3. Merge sorted chunks into Zone 3.
        Self::run_from_sorted_chunks(config, &header, sorted_chunks)
    }

    /// Run Zone 3 from pre-sorted chunks (sort-accumulate already done).
    ///
    /// Spawns the K-way merge on a feeder thread and immediately joins the calling thread
    /// to the Zone 3 worker pool. This is the shared tail of both `run_from_sort` (which
    /// reads from a BAM file) and direct callers that sort-accumulate from an iterator.
    ///
    /// # Errors
    ///
    /// Returns an error if any stage fails, I/O fails, or the scheduler encounters an error.
    pub fn run_from_sorted_chunks(
        config: &PipelineConfig,
        header: &noodles::sam::Header,
        sorted_chunks: fgumi_lib::sort::SortedChunks,
    ) -> Result<()> {
        // 1. Build Zone 3 queues and stages, start scheduler (including write stage).
        let zone3 = Zone3::start(config, header)?;

        // 2. Spawn merge on a dedicated feeder thread so the calling thread can join
        //    the Zone 3 worker pool immediately.
        let q_position_batch = zone3.q_position_batch.clone();
        let merge_handle = thread::spawn(move || -> Result<u64> {
            // Guard: close the queue when this scope exits (success, error, or panic)
            // to ensure the close cascade always fires and wait_for_completion won't hang.
            let close_guard = CloseOnDrop(q_position_batch.clone());
            let mut batcher = PositionBatcher::new(q_position_batch);
            let records_merged =
                merge_sorted_chunks_to_sink(sorted_chunks, &mut batcher).context("merge phase")?;
            batcher.finish().context("flushing final position batch")?;
            let total_batches = batcher.batches_pushed();
            info!("Pipeline: merge complete, fed {total_batches} position batches to Zone 3");
            drop(close_guard); // explicit close for clarity
            Ok(records_merged)
        });

        // 3. Calling thread immediately joins the pool as a worker.
        //    The merge thread will close q_position_batch (via batcher drop / close cascade)
        //    when done, which cascades through all queues to signal completion.
        zone3.wait_for_completion(0)?;

        // 4. Wait for merge thread and propagate any errors.
        let _records_merged =
            merge_handle.join().map_err(|e| anyhow::anyhow!("merge thread panicked: {e:?}"))??;

        info!("Pipeline: done");
        Ok(())
    }

    /// Run the full two-phase pipeline: Phase A (extract → align → zipper → sort)
    /// followed by Phase B (Zone 3: group → consensus → filter → compress → write).
    ///
    /// Phase A uses a work-stealing pool for Extract, `ToFastq`, and Zipper stages.
    /// Two dedicated I/O threads handle the aligner stdout (SAM parsing) and sort
    /// accumulation (draining `q_zippered`). These I/O threads do not count against
    /// `--threads`.
    ///
    /// The close cascade propagates through the queues:
    /// - Extract done → close `q_extracted`
    /// - `q_extracted` drained by `ToFastq` → close aligner stdin + close `q_unmapped`
    /// - Aligner finishes → stdout EOF → stdout-reader closes `q_mapped`
    /// - Zipper thread drains `q_unmapped` + `q_mapped` → pushes to `q_zippered`
    /// - After Zipper thread exits → close `q_zippered`
    /// - `q_zippered` drained → sort-accumulate produces `SortedChunks`
    /// - `SortedChunks` merged into Zone 3
    ///
    /// # Errors
    ///
    /// Returns an error if any phase, stage, or I/O thread fails.
    #[expect(clippy::too_many_lines, reason = "orchestration method wires multiple phases")]
    pub fn run_two_phase(config: TwoPhaseConfig) -> Result<()> {
        info!("Pipeline: starting two-phase pipeline (Phase A + Phase B)");

        // 1. Build Phase A graph (output_header is not yet available — it will be set
        //    by the stdout-reader thread after reading the SAM header from the aligner).
        let graph = Arc::new(PhaseAGraph::new(
            config.extract_source,
            config.extract_params,
            config.quality_encoding,
            Arc::clone(&config.tag_info),
            config.skip_pa_tags,
            config.merge_fn,
            config.batch_state,
            config.aligner_stdin,
            Arc::clone(&config.cancel),
            config.correction_config,
        ));

        // 2. Spawn stdout-reader thread (dedicated, not counted in --threads).
        //    Reads the SAM header from aligner stdout, builds the output header,
        //    sets it on the graph, sends it to the main thread, then reads records.
        let (header_tx, header_rx) = std::sync::mpsc::sync_channel::<Arc<noodles::sam::Header>>(1);
        let q_mapped = graph.q_mapped.clone();
        let graph_for_stdout = Arc::clone(&graph);
        let build_output_header_fn = config.build_output_header_fn;
        let aligner_stdout = config.aligner_stdout;
        let stdout_cancel = Arc::clone(&config.cancel);
        let stdout_reader = thread::Builder::new()
            .name("stdout-reader".into())
            .spawn(move || -> Result<u64> {
                // Read the SAM header from aligner stdout.
                let buf_reader = std::io::BufReader::with_capacity(256 * 1024, aligner_stdout);
                let mut sam_reader = noodles::sam::io::Reader::new(buf_reader);
                let mapped_header = sam_reader.read_header()?;
                info!("stdout-reader: read SAM header from aligner");

                // Build the final output header and publish it.
                let output_header = build_output_header_fn(&mapped_header)?;
                let output_header = Arc::new(output_header);
                let _ = graph_for_stdout.output_header.set(Arc::clone(&output_header));
                let _ = header_tx.send(Arc::clone(&output_header));

                // Read records from aligner stdout.
                let iter = TemplateIterator::new(
                    sam_reader.record_bufs(&mapped_header).map(|r| r.map_err(anyhow::Error::from)),
                );
                let mut count = 0u64;
                for template_result in iter {
                    if stdout_cancel.load(Ordering::Acquire) {
                        break;
                    }
                    let template = template_result?;
                    let mem = template.estimate_heap_size();
                    q_mapped.push_blocking(template, mem)?;
                    count += 1;
                }
                q_mapped.close();
                info!("stdout-reader: finished, pushed {count} templates to q_mapped");
                Ok(count)
            })
            .context("failed to spawn stdout-reader thread")?;

        // 3. Spawn dedicated Zipper thread.
        //    Sequentially matches q_unmapped + q_mapped templates by queryname
        //    (guaranteed to be in the same order), merges, encodes to BAM bytes,
        //    and pushes to q_zippered. No HashMap needed — ordering is guaranteed
        //    because ToFastq is exclusive and bwa preserves input order.
        let graph_for_zipper = Arc::clone(&graph);
        let zipper_handle = thread::Builder::new()
            .name("zipper".into())
            .spawn(move || -> Result<u64> { phase_a_steps::zipper_thread_loop(&graph_for_zipper) })
            .context("failed to spawn zipper thread")?;

        // 4. Spawn sort-accumulate thread (dedicated, not counted in --threads).
        //    Drains q_zippered and accumulates records into SortedChunks.
        let q_zippered = graph.q_zippered.clone();
        let sort_memory_limit = config.sort_memory_limit;
        let sort_temp_dir = config.sort_temp_dir.clone();
        let sort_threads = config.sort_threads;
        let graph_for_sort = Arc::clone(&graph);
        let sort_cancel = Arc::clone(&config.cancel);
        let sort_handle = thread::Builder::new()
            .name("sort-accumulate".into())
            .spawn(move || -> Result<fgumi_lib::sort::SortedChunks> {
                // Wait for the output header to become available (set by stdout-reader).
                let sort_header = loop {
                    if let Some(h) = graph_for_sort.output_header.get() {
                        break h.as_ref().clone();
                    }
                    if sort_cancel.load(Ordering::Acquire) {
                        anyhow::bail!("sort-accumulate cancelled while waiting for header");
                    }
                    thread::sleep(std::time::Duration::from_millis(1));
                };

                // Build an iterator that drains q_zippered via blocking pop.
                let iter = std::iter::from_fn(|| {
                    if sort_cancel.load(Ordering::Acquire) {
                        return None;
                    }
                    let item = q_zippered.pop()?;
                    let mem = item.len();
                    q_zippered.complete_in_flight();
                    q_zippered.release_memory(mem);
                    Some(item)
                });

                let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                    .memory_limit(sort_memory_limit)
                    .threads(sort_threads)
                    .temp_compression(1);
                if let Some(ref temp_dir) = sort_temp_dir {
                    sorter = sorter.temp_dir(temp_dir.clone());
                }
                let sorted_chunks = sorter
                    .sort_accumulate_from_iter(iter, &sort_header)
                    .context("sort-accumulate from q_zippered")?;
                info!(
                    "sort-accumulate: complete ({} records, {} disk chunks)",
                    sorted_chunks.stats.total_records, sorted_chunks.stats.chunks_written
                );
                Ok(sorted_chunks)
            })
            .context("failed to spawn sort-accumulate thread")?;

        // 5. Run Phase A worker pool (Extract + optional Correct + ToFastq).
        //    Zipper and SortAccumulate are dedicated threads.
        let mut active_steps = [false; PhaseAStep::NUM_STEPS];
        active_steps[PhaseAStep::Extract as usize] = true;
        active_steps[PhaseAStep::Correct as usize] = graph.correction_config.is_some();
        active_steps[PhaseAStep::ToFastq as usize] = true;

        let num_workers = config.threads.max(1);
        let umi_cache_size = graph.correction_config.as_ref().map_or(0, |c| c.cache_size);
        let (handles, worker_slots, error_rx) =
            start_phase_a_pool(&graph, num_workers, active_steps, &config.cancel, umi_cache_size);
        info!("Pipeline Phase A: started {num_workers} worker threads");

        // Wait for all Phase A workers to finish.
        wait_phase_a_pool(handles, &worker_slots, &error_rx)
            .context("Phase A worker pool error")?;
        info!("Pipeline Phase A: all workers finished");

        // Close q_unmapped and aligner stdin. The pool may have exited before
        // ToFastq got a chance to close these (the last ToFastq step processes
        // an item and returns Success; the pool then exits via is_pool_done()
        // before ToFastq runs again to see the empty/closed q_extracted).
        graph.close_aligner_stdin();
        graph.q_unmapped.close();
        if let Some(ref q) = graph.q_corrected {
            q.close();
        }

        // 6. Wait for aligner process to complete.
        config.aligner_process.wait().context("aligner subprocess failed")?;
        info!("Pipeline Phase A: aligner process exited successfully");

        // 7. Wait for stdout-reader thread.
        let mapped_count = stdout_reader
            .join()
            .map_err(|e| anyhow::anyhow!("stdout-reader thread panicked: {e:?}"))??;
        info!("Pipeline Phase A: stdout-reader pushed {mapped_count} templates");

        // 8. Wait for Zipper thread (closes q_zippered when done).
        let zipper_count = zipper_handle
            .join()
            .map_err(|e| anyhow::anyhow!("zipper thread panicked: {e:?}"))??;
        info!("Pipeline Phase A: zipper merged {zipper_count} records");

        // Close q_zippered to signal sort-accumulate to finish.
        graph.q_zippered.close();

        // 9. Wait for sort-accumulate thread.
        let sorted_chunks = sort_handle
            .join()
            .map_err(|e| anyhow::anyhow!("sort-accumulate thread panicked: {e:?}"))??;
        info!(
            "Pipeline Phase A: sort-accumulate produced {} records in {} chunks",
            sorted_chunks.stats.total_records, sorted_chunks.stats.chunks_written
        );

        // 10. Phase B: merge sorted chunks through Zone 3.
        info!("Pipeline Phase B: starting Zone 3 from sorted chunks");
        let output_header = header_rx
            .recv()
            .map_err(|_| anyhow::anyhow!("failed to receive output header from stdout-reader"))?;
        let phase_b_header = output_header.as_ref().clone();
        Self::run_from_sorted_chunks(&config.pipeline_config, &phase_b_header, sorted_chunks)?;

        info!("Pipeline: two-phase pipeline complete");
        Ok(())
    }
}

// ============================================================================
// Zone3 — internal orchestration of the Zone 3 stages
// ============================================================================

/// Holds the running Zone 3 pipeline: queues, scheduler, and shared state.
///
/// All threads participate in work-stealing, including the write stage. There is no
/// dedicated writer thread; writing is the highest-priority stage for the write-preferred
/// worker.
struct Zone3 {
    /// Queue where position batches enter Zone 3.
    q_position_batch: WorkQueue<SequencedItem<PositionGroupBatch>>,
    /// Final queue — closed+empty signals all work has cascaded through the pipeline.
    q_compressed: WorkQueue<SequencedItem<CompressedBatch>>,
    /// The scheduler managing worker threads.
    scheduler: Scheduler,
    /// Worker thread handles.
    worker_handles: Vec<thread::JoinHandle<()>>,
    /// Per-worker state (for stats collection after join).
    worker_states: Vec<WorkerStateSlot>,
    /// Number of workers (for creating the calling thread's `WorkerState`).
    num_workers: usize,
    /// Shared write state (file handle, reorder buffer, batch counter).
    write_state: Arc<Mutex<WriteState>>,
    /// Per-stage timing stats.
    stage_stats: Vec<Arc<StageStats>>,
}

impl Zone3 {
    /// Create queues, stages, and start the scheduler with the typed `StageGraph`.
    ///
    /// The output file is created and the BAM header is written before workers start,
    /// so that the write runner only needs to append pre-compressed BGZF blocks.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created or the header cannot be written.
    fn start(config: &PipelineConfig, header: &noodles::sam::Header) -> Result<Self> {
        let cancel = Arc::new(AtomicBool::new(false));

        // -- Queues --
        let q_position_batch: WorkQueue<SequencedItem<PositionGroupBatch>> =
            WorkQueue::new(QUEUE_CAPACITY, QUEUE_MEMORY_LIMIT, "position_batch");
        let q_mi_batch: WorkQueue<SequencedItem<Vec<MiGroup>>> =
            WorkQueue::new(QUEUE_CAPACITY, QUEUE_MEMORY_LIMIT, "mi_batch");
        let q_consensus: WorkQueue<SequencedItem<fgumi_consensus::caller::ConsensusOutput>> =
            WorkQueue::new(QUEUE_CAPACITY, QUEUE_MEMORY_LIMIT, "consensus");
        let q_filtered: WorkQueue<SequencedItem<Vec<u8>>> =
            WorkQueue::new(QUEUE_CAPACITY, QUEUE_MEMORY_LIMIT, "filtered");
        let q_compressed: WorkQueue<SequencedItem<CompressedBatch>> =
            WorkQueue::new(QUEUE_CAPACITY, QUEUE_MEMORY_LIMIT, "compressed");

        // -- Open output file and write BAM header before starting workers --
        let mut file = std::fs::File::create(&config.output)
            .with_context(|| format!("creating output BAM: {}", config.output.display()))?;
        write_bam_header_to_file(&mut file, header, config.compression_level)?;

        let write_state = Arc::new(Mutex::new(WriteState {
            file,
            reorder_buf: ReorderBuffer::new(),
            batches_written: 0,
        }));

        // -- Stages --
        let group_assign_stage = Arc::new(GroupAndAssignStage::new(
            Arc::clone(&config.assigner),
            config.umi_tag,
            config.assign_tag,
            config.group_filter_config.clone(),
        ));
        let consensus_stage = Arc::new(ConsensusStage::new(Arc::clone(&config.caller_factory)));
        let filter_stage = Arc::new(FilterStage::new(config.filter_config.clone()));
        let compress_stage = Arc::new(CompressStage::new(config.compression_level));

        // -- Per-stage stats --
        let group_assign_stats = Arc::new(StageStats::new("group_assign"));
        let consensus_stats = Arc::new(StageStats::new("consensus"));
        let filter_stats = Arc::new(StageStats::new("filter"));
        let compress_stats = Arc::new(StageStats::new("compress"));
        let writer_stats = Arc::new(StageStats::new("write"));
        let stage_stats: Vec<Arc<StageStats>> = vec![
            Arc::clone(&group_assign_stats),
            Arc::clone(&consensus_stats),
            Arc::clone(&filter_stats),
            Arc::clone(&compress_stats),
            Arc::clone(&writer_stats),
        ];

        // -- Build StageGraph --
        let graph = StageGraph {
            group_assign: group_assign_stage,
            consensus: consensus_stage,
            filter: filter_stage,
            compress: compress_stage,
            q_position_batch: q_position_batch.clone(),
            q_mi_batch,
            q_consensus,
            q_filtered,
            q_compressed: q_compressed.clone(),
            write_state: Arc::clone(&write_state),
            stats: [
                group_assign_stats,
                consensus_stats,
                filter_stats,
                compress_stats,
                writer_stats,
            ],
            cancel,
        };

        // -- Scheduler --
        let num_workers = config.threads.max(1);
        let scheduler = Scheduler::new(graph);
        let (worker_handles, worker_states) = scheduler.start(num_workers);
        info!("Pipeline: started {num_workers} worker threads (unified model)");

        Ok(Self {
            q_position_batch,
            q_compressed,
            scheduler,
            worker_handles,
            worker_states,
            num_workers,
            write_state,
            stage_stats,
        })
    }

    /// Wait for all input to be fully processed through the pipeline.
    ///
    /// The calling thread participates as a worker, running the same priority
    /// loop as spawned workers. This avoids idle spinning and ensures progress even
    /// with `--threads 1`.
    ///
    /// **Important:** The caller must ensure that `q_position_batch` will be closed
    /// (either before calling this method, or from a feeder thread) so that the
    /// close cascade can propagate through all queues and signal completion.
    #[expect(clippy::too_many_lines, reason = "inline worker loop mirrors spawned worker logic")]
    fn wait_for_completion(self, _expected_position_batches: u64) -> Result<()> {
        // The calling thread participates as an additional worker.
        let graph = Arc::clone(self.scheduler.graph());
        let _cancel_flag = self.scheduler.cancel_flag();
        let (_caller_error_tx, _caller_error_rx) = crossbeam_channel::unbounded::<anyhow::Error>();
        let caller_worker_id = self.num_workers; // one past the last spawned worker
        let total_workers = self.num_workers + 1;
        let mut caller_worker = WorkerState::new(caller_worker_id, total_workers);

        loop {
            if self.scheduler.is_cancelled() {
                break;
            }
            if let Some(err) = self.scheduler.check_error() {
                self.scheduler.cancel();
                return Err(err).context("scheduler error during pipeline execution");
            }

            // Run one iteration of the worker loop inline.
            let mut did_work = false;

            // 1. Try to push held items.
            did_work |=
                crate::pipeline::scheduler::try_push_held_compress_pub(&graph, &mut caller_worker);
            did_work |=
                crate::pipeline::scheduler::try_push_held_filter_pub(&graph, &mut caller_worker);
            did_work |=
                crate::pipeline::scheduler::try_push_held_consensus_pub(&graph, &mut caller_worker);
            did_work |= crate::pipeline::scheduler::try_push_held_group_assign_pub(
                &graph,
                &mut caller_worker,
            );

            // 2. Sample backpressure and get priorities.
            let bp = BackpressureState::sample(&graph);
            let priorities = *caller_worker.scheduler.get_priorities(&bp);

            // 3. Execute one stage in priority order.
            for &stage in &priorities {
                if self.scheduler.is_cancelled() {
                    break;
                }

                let result = match stage {
                    PipelineStep::Write => {
                        crate::pipeline::scheduler::try_step_write_pub(&graph, &mut caller_worker)
                    }
                    PipelineStep::Compress => crate::pipeline::scheduler::try_step_compress_pub(
                        &graph,
                        &mut caller_worker,
                    ),
                    PipelineStep::Filter => {
                        crate::pipeline::scheduler::try_step_filter_pub(&graph, &mut caller_worker)
                    }
                    PipelineStep::Consensus => crate::pipeline::scheduler::try_step_consensus_pub(
                        &graph,
                        &mut caller_worker,
                    ),
                    PipelineStep::GroupAssign => {
                        crate::pipeline::scheduler::try_step_group_assign_pub(
                            &graph,
                            &mut caller_worker,
                        )
                    }
                };

                match result {
                    StepResult::Success => {
                        caller_worker.scheduler.record_outcome(stage, true);
                        did_work = true;
                        break;
                    }
                    StepResult::OutputFull => {
                        caller_worker.scheduler.record_outcome(stage, false);
                    }
                    StepResult::InputEmpty => {}
                    StepResult::Error(e) => {
                        self.scheduler.cancel();
                        return Err(e).context("error in calling thread worker loop");
                    }
                }
            }

            if did_work {
                caller_worker.reset_backoff();
            } else {
                // Check completion.
                let done = self.q_compressed.is_closed_and_empty()
                    && self.write_state.lock().expect("lock").reorder_buf.is_empty();

                if done {
                    break;
                }

                caller_worker.sleep_backoff();
                caller_worker.increase_backoff();
            }
        }

        // Signal scheduler to stop and wait for workers.
        self.scheduler.cancel();
        self.scheduler.wait(self.worker_handles, &self.worker_states)?;

        // Log caller worker stats.
        log::info!("  Caller thread:");
        caller_worker.log_stats();

        // Write BGZF EOF marker and flush.
        let mut state = self.write_state.lock().expect("write state lock poisoned");
        state.file.write_all(&BGZF_EOF).context("writing BGZF EOF")?;
        state.file.flush().context("flushing output BAM")?;

        info!("Pipeline writer: finished, wrote {} batches", state.batches_written);

        // Log per-stage timing stats.
        info!("Pipeline stage stats:");
        #[expect(clippy::cast_precision_loss, reason = "item counts are approximate for display")]
        for stats in &self.stage_stats {
            let secs = stats.total_secs();
            let items = stats.items();
            let rate = if secs > 0.0 { items as f64 / secs } else { 0.0 };
            info!(
                "  {:<14} {:.2}s  {:>10} items  {:>10.0} items/sec",
                stats.name(),
                secs,
                items,
                rate
            );
        }

        Ok(())
    }
}

// ============================================================================
// CloseOnDrop — ensures a WorkQueue is closed even on panic/error
// ============================================================================

/// RAII guard that closes a [`WorkQueue`] when dropped.
///
/// Used by the merge feeder thread to guarantee the close cascade fires,
/// even if the merge errors or the thread panics.
struct CloseOnDrop<T>(WorkQueue<SequencedItem<T>>);

impl<T> Drop for CloseOnDrop<T> {
    fn drop(&mut self) {
        self.0.close();
    }
}

// ============================================================================
// BAM header writing
// ============================================================================

/// Write the BAM header to a file using BGZF compression.
///
/// Uses a temporary in-memory BGZF compressor to produce BGZF-compressed header
/// blocks, then writes those blocks directly to the file.
fn write_bam_header_to_file(
    file: &mut std::fs::File,
    header: &noodles::sam::Header,
    compression_level: u32,
) -> Result<()> {
    use fgumi_bgzf::writer::InlineBgzfCompressor;

    let mut compressor = InlineBgzfCompressor::new(compression_level);

    // Serialize the BAM header into the compressor.
    // BAM magic
    compressor.write_all(fgumi_raw_bam::BAM_MAGIC)?;

    // Header text
    let mut sam_writer = noodles::sam::io::Writer::new(Vec::new());
    sam_writer.write_header(header)?;
    let header_bytes = sam_writer.into_inner();
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "BAM header text is always < 2 GB"
    )]
    let l_text = header_bytes.len() as i32;
    compressor.write_all(&l_text.to_le_bytes())?;
    compressor.write_all(&header_bytes)?;

    // Reference sequences
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "number of reference sequences always fits in i32"
    )]
    let n_ref = header.reference_sequences().len() as i32;
    compressor.write_all(&n_ref.to_le_bytes())?;

    for (name, map) in header.reference_sequences() {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "reference name length always fits in u32"
        )]
        let l_name = (name.len() + 1) as u32;
        compressor.write_all(&l_name.to_le_bytes())?;
        compressor.write_all(name.as_ref())?;
        compressor.write_all(&[0u8])?; // NUL terminator

        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_possible_wrap,
            reason = "reference length always fits in i32"
        )]
        let l_ref = map.length().get() as i32;
        compressor.write_all(&l_ref.to_le_bytes())?;
    }

    compressor.flush()?;
    let blocks = compressor.take_blocks();
    for block in &blocks {
        file.write_all(&block.data).context("writing header BGZF block")?;
    }

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stage_enum_variants() {
        assert_eq!(Stage::Sort, Stage::Sort);
        assert_eq!(Stage::Group, Stage::Group);
        assert_ne!(Stage::Sort, Stage::Group);
    }

    #[test]
    fn test_queue_defaults_are_reasonable() {
        const { assert!(QUEUE_CAPACITY > 0) };
        const { assert!(QUEUE_MEMORY_LIMIT > 0) };
    }
}

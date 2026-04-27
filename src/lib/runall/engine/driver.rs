//! # Pipeline Driver
//!
//! Orchestrates pipeline execution: spawns the work-stealing pool,
//! threads [`SpecialStage`](super::special_stage::SpecialStage)
//! instances onto dedicated OS threads, collects stats, drives
//! metrics output, and handles cancellation.
//!
//! ## Topology
//!
//! ```text
//! Source -> [Stage | Stage | SpecialStage | Stage | ...] -> Sink
//!            \-------- work-stealing pool --------/
//!            special stages on dedicated OS threads
//! ```
//!
//! Each `Stage -> Stage` transition is a
//! [`StageQueue<T>`](super::queue::StageQueue) backed by
//! `crossbeam_queue::ArrayQueue`. Items flow between stages as
//! [`SequencedItem<T>`](super::stage::SequencedItem) values carrying a
//! monotonic `ordinal` and a memory estimate so downstream reassembly
//! (e.g. [`Coalesce`](super::stages::coalesce::Coalesce),
//! [`SortStage`](super::stages::sort::SortStage), and the
//! [`ReorderBuffer`](super::reorder::ReorderBuffer) inside
//! `BamFileWrite`) can restore order.
//!
//! ## Threading
//!
//! - **Pool workers.** N workers (configured by
//!   [`PipelineConfig::worker_threads`]) cooperatively service every
//!   non-special stage via the work-stealing [`Pool`]. Workers are
//!   spawned through [`run_worker`].
//! - **Special stages.** Each
//!   [`SpecialStage`](super::special_stage::SpecialStage) owns a
//!   dedicated OS thread; it drains its input queue and drives its
//!   output queue directly, outside the pool. Use `SpecialStage`
//!   when a stage needs single-threaded ordering
//!   (`PositionBatchStage`, `Coalesce`), barrier semantics
//!   (`SortStage`), or its own inner thread coordination
//!   (`AlignAndMerge`).
//! - **Panic catch.** Every worker body is wrapped in `catch_unwind`
//!   via `run_and_catch_panic_named` so panics convert to
//!   `anyhow::Error` with stage-name context instead of hanging
//!   siblings on a queue that never drains.
//!
//! ## Cancellation
//!
//! A shared [`CancelToken`] propagates SIGINT (and pipeline-internal
//! errors) across all threads. Workers
//! and special stages check the token inside push/pop loops; once
//! triggered, stages finish the record they hold and exit. Sinks that
//! use an `AtomicOutputGuard` discard their `.tmp` output rather than
//! promoting it to the final path on cancel.
//!
//! ## Memory gating
//!
//! A shared [`MemoryTracker`] enforces a global byte ceiling on
//! queued batches, using the per-item
//! `mem_bytes` estimate each stage supplies. See
//! [`crate::runall::engine::memory`] for scope and accounting rules.
//!
//! ## Deadlock detection
//!
//! A watchdog ([`watch`]) periodically samples per-stage counters
//! ([`CountingInputQueue`] / [`CountingOutputQueue`]) and cancels the
//! pipeline if no stage makes progress within the configured window.

use std::any::Any;
use std::marker::PhantomData;
use std::path::Path;
use std::sync::Arc;
use std::thread;

use anyhow::{Context, Result};

use std::sync::atomic::AtomicU64;

use super::cancel::CancelToken;
use super::counting::{CountingInputQueue, CountingOutputQueue};
use super::deadlock::{QueueProgress, WatchdogConfig, WatchdogOutcome, watch};
use super::erased_queue::{ErasingOutputQueue, UnerasingInputQueue};
use super::memory::MemoryTracker;
use super::pool::{Pool, StageInvoker};
use super::queue::StageQueue;
use super::sink::{InputQueue, Sink};
use super::source::{OutputQueue, Source};
use super::special_stage::ErasedSpecialStage;
use super::stage::Stage;
use super::stage_source::ErasedStageSource;
use super::stats::PipelineStats;
use super::test_support;
use super::worker::run_worker;

/// Extract a human-readable message from a panic payload.
///
/// Rust's `thread::Result::Err(Box<dyn Any + Send>)` carries the panic
/// value. For `panic!("literal")` the payload is `&'static str`; for
/// `panic!("{}", value)` it is `String`. Using `{:?}` on the `Box<dyn Any>`
/// stringifies the box itself and hides the actual message. Downcasting
/// to those common shapes recovers the original string; unknown types
/// fall back to a `TypeId` descriptor so at least the type is visible.
pub(crate) fn extract_panic_message(payload: &(dyn Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&'static str>() {
        return (*s).to_string();
    }
    if let Some(s) = payload.downcast_ref::<String>() {
        return s.clone();
    }
    format!("(non-string panic payload: {:?})", payload.type_id())
}

/// Join a thread handle, propagating the first error observed and cancelling
/// the pipeline on failure so siblings can exit their retry loops.
///
/// If `handle` panicked, siblings would otherwise block forever on a queue
/// that never drains; cancelling signals them to exit promptly. We still
/// return the first error seen so the caller can propagate it.
///
/// Panic payloads are unwrapped via [`extract_panic_message`] so the
/// caller-visible error contains the actual `panic!()` message instead of
/// a `Box<dyn Any>` debug stringification.
fn join_thread_or_cancel<T>(
    handle: thread::JoinHandle<Result<T>>,
    cancel: &CancelToken,
) -> Result<T> {
    let name = handle.thread().name().unwrap_or("unnamed").to_string();
    match handle.join() {
        Ok(Ok(value)) => Ok(value),
        Ok(Err(e)) => {
            cancel.cancel();
            Err(e)
        }
        Err(payload) => {
            cancel.cancel();
            let msg = extract_panic_message(payload.as_ref());
            Err(anyhow::anyhow!("thread '{name}' panicked: {msg}"))
        }
    }
}

/// Run `body` while catching any panic. On panic, mark the shared cancel
/// token so sibling threads exit their retry loops; then translate the
/// panic payload into an `Err(anyhow::Error)` so `join_thread_or_cancel`
/// reports it as an ordinary error instead of waiting for every other
/// thread to finish first.
///
/// Workers use this because a panic deep inside `Stage::process` would
/// otherwise leave siblings (source, sink, neighbouring stages) blocked
/// forever on a queue that never drains.
///
/// We wrap `body` in [`std::panic::AssertUnwindSafe`] because worker
/// bodies hold `Box<dyn Stage>` / `Box<dyn ErasedStage>` which are not
/// `UnwindSafe`. That is sound here: after a worker panics, the cancel
/// token is set and no shared state is re-accessed by siblings except
/// through the `StageQueue` (which is inherently `Sync`).
pub(crate) fn run_and_catch_panic<F>(cancel: &CancelToken, body: F) -> Result<()>
where
    F: FnOnce() -> Result<()>,
{
    run_and_catch_panic_named(cancel, None, body)
}

/// Variant of [`run_and_catch_panic`] that attributes a caught panic to a
/// named pipeline stage. Use this when spawning a thread that owns a single
/// named stage (e.g. a special/self-threaded stage); the stage name is
/// included in the synthesized error message so operators can tell which
/// stage blew up without reading thread names out of log lines.
///
/// When `stage_name` is `None`, the error reads identically to
/// `run_and_catch_panic` for backwards compatibility.
pub(crate) fn run_and_catch_panic_named<F>(
    cancel: &CancelToken,
    stage_name: Option<&str>,
    body: F,
) -> Result<()>
where
    F: FnOnce() -> Result<()>,
{
    let wrapped = std::panic::AssertUnwindSafe(body);
    match std::panic::catch_unwind(wrapped) {
        Ok(Ok(())) => Ok(()),
        Ok(Err(e)) => {
            cancel.cancel();
            Err(e)
        }
        Err(payload) => {
            cancel.cancel();
            let msg = extract_panic_message(payload.as_ref());
            match stage_name {
                Some(name) => Err(anyhow::anyhow!("stage '{name}' panicked: {msg}")),
                None => Err(anyhow::anyhow!("worker thread panicked: {msg}")),
            }
        }
    }
}

/// Trait object view of a queue that can be closed.
///
/// Used by `QueueCloseGuard` so a heterogeneous set of `StageQueue<T>` with
/// differing `T` can share a single close-on-drop list.
trait Closeable: Send + Sync {
    fn close_queue(&self);
}

impl<T: Send> Closeable for StageQueue<T> {
    fn close_queue(&self) {
        self.close();
    }
}

/// Scope guard that closes every registered queue when dropped, unless
/// `disarmed` has been called.
///
/// Purpose: if the driver unwinds (panic, early-return) before manually
/// closing every intermediate queue, any thread still spinning on a retry
/// loop would block forever. The guard ensures queues are closed on every
/// exit path, including panics during join/close.
struct QueueCloseGuard {
    queues: Vec<Arc<dyn Closeable>>,
    disarmed: bool,
}

impl QueueCloseGuard {
    fn new() -> Self {
        Self { queues: Vec::new(), disarmed: false }
    }

    fn register<T: Send + 'static>(&mut self, queue: Arc<StageQueue<T>>) {
        let as_closeable: Arc<dyn Closeable> = queue;
        self.queues.push(as_closeable);
    }

    fn disarm(&mut self) {
        self.disarmed = true;
    }
}

impl Drop for QueueCloseGuard {
    fn drop(&mut self) {
        if self.disarmed {
            return;
        }
        for q in &self.queues {
            q.close_queue();
        }
    }
}

/// Minimal configuration for the pipeline driver.
///
/// Defaults are tuned for BAM-scale workloads: 2 GB per-queue and 4 GB global
/// memory budgets, matching the v1 `unified_pipeline` tuning that was
/// settled after iterated benchmarking. Small-workload tests can override
/// these to keep footprints tight.
#[derive(Debug, Clone, Copy)]
pub struct PipelineConfig {
    /// Number of worker threads for stages.
    pub worker_threads: usize,
    /// Capacity (slots) per inter-stage queue.
    pub queue_capacity: usize,
    /// Memory limit (bytes) per inter-stage queue.
    pub queue_memory_limit: usize,
    /// Global memory limit (bytes) across all queues.
    pub global_memory_limit: usize,
}

/// Default per-queue memory limit (bytes): 2 GiB.
pub const DEFAULT_QUEUE_MEMORY_LIMIT: usize = 2 * 1024 * 1024 * 1024;

/// Default global memory limit (bytes): 4 GiB.
pub const DEFAULT_GLOBAL_MEMORY_LIMIT: usize = 4 * 1024 * 1024 * 1024;

/// Default queue capacity (slot count).
pub const DEFAULT_QUEUE_CAPACITY: usize = 1024;

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            worker_threads: 1,
            queue_capacity: DEFAULT_QUEUE_CAPACITY,
            queue_memory_limit: DEFAULT_QUEUE_MEMORY_LIMIT,
            global_memory_limit: DEFAULT_GLOBAL_MEMORY_LIMIT,
        }
    }
}

/// Lower bound on auto-tuned queue capacity (slots).
pub const AUTO_TUNED_QUEUE_CAPACITY_MIN: usize = 64;

/// Upper bound on auto-tuned queue capacity (slots).
pub const AUTO_TUNED_QUEUE_CAPACITY_MAX: usize = 256;

/// Slot-per-thread multiplier used by [`PipelineConfig::auto_tuned`].
///
/// Chosen so that even a single-threaded configuration starts at the lower
/// bound (1 * 16 = 16 clamps up to 64), and 16 threads saturates the upper
/// bound (16 * 16 = 256). Matches the v1 `unified_pipeline` BAM-pipeline
/// tuning.
pub const AUTO_TUNED_SLOTS_PER_THREAD: usize = 16;

impl PipelineConfig {
    /// Build a config with queue capacity scaled to thread count.
    ///
    /// `queue_capacity = clamp(threads * 16, 64, 256)`. Other fields use
    /// the same defaults as [`PipelineConfig::default`]. Recommended entry
    /// point when the caller already knows its worker-thread budget; raw
    /// `default()` is retained so test code and ad-hoc callers don't have
    /// to pass `threads` explicitly.
    ///
    /// Matches the v1 BAM-pipeline tuning: threads<=4 get 64 slots (the
    /// lower bound), threads>=16 get 256 slots (the upper bound), and
    /// values in between scale linearly.
    #[must_use]
    pub fn auto_tuned(threads: usize) -> Self {
        let capacity = (threads * AUTO_TUNED_SLOTS_PER_THREAD)
            .clamp(AUTO_TUNED_QUEUE_CAPACITY_MIN, AUTO_TUNED_QUEUE_CAPACITY_MAX);
        Self {
            worker_threads: threads,
            queue_capacity: capacity,
            queue_memory_limit: DEFAULT_QUEUE_MEMORY_LIMIT,
            global_memory_limit: DEFAULT_GLOBAL_MEMORY_LIMIT,
        }
    }
}

/// Run a single-stage pipeline: Source -> Stage -> Sink.
///
/// Spawns one source thread, `worker_threads` worker threads, and one sink thread.
/// Blocks until all threads complete.
///
/// # Errors
///
/// Returns an error if any thread fails or panics.
pub fn run_single_stage<S, St, Sk>(
    source: Box<S>,
    stage_factory: impl Fn() -> St + Send + Sync + 'static,
    sink: Box<Sk>,
    cancel: &CancelToken,
    cfg: PipelineConfig,
) -> Result<()>
where
    S: Source + ?Sized + 'static,
    St: Stage<Input = S::Output, Output = Sk::Input> + 'static,
    Sk: Sink + ?Sized + 'static,
{
    let global_memory = Arc::new(MemoryTracker::new(cfg.global_memory_limit));

    let q_in: Arc<StageQueue<S::Output>> = Arc::new(StageQueue::new(
        "source_out",
        cfg.queue_capacity,
        cfg.queue_memory_limit,
        global_memory.clone(),
    ));
    let q_out: Arc<StageQueue<Sk::Input>> = Arc::new(StageQueue::new(
        "stage_out",
        cfg.queue_capacity,
        cfg.queue_memory_limit,
        global_memory.clone(),
    ));

    // On any error path — including panic during driver logic — ensure every
    // intermediate queue is closed so threads retrying pushes can exit. We
    // disarm on the happy path just before returning Ok.
    let mut close_guard = QueueCloseGuard::new();
    close_guard.register(q_in.clone());
    close_guard.register(q_out.clone());

    // Counters for source-produced and sink-consumed items. After all threads
    // join we assert `produced == consumed` to catch silent data loss in a
    // stage, scheduler, or panic handler.
    let produced = Arc::new(AtomicU64::new(0));
    let consumed = Arc::new(AtomicU64::new(0));

    // Source thread. Wraps the typed queue in a `CountingOutputQueue` so
    // every successful push increments the produced counter. No boxing of
    // items occurs on the hot push path (`push` stays typed).
    let q_in_for_source = q_in.clone();
    let cancel_for_source = cancel.clone();
    let source_out: Box<dyn OutputQueue<S::Output>> = Box::new(CountingOutputQueue::new(
        Box::new(q_in_for_source) as Box<dyn OutputQueue<S::Output>>,
        produced.clone(),
    ));
    let source_handle = thread::Builder::new()
        .name("pipeline_source".into())
        .spawn(move || source.run(source_out, cancel_for_source))
        .context("failed to spawn source thread")?;

    // Worker threads. Wrap each worker body in `catch_unwind` so a panic
    // inside Stage::process cancels the pipeline and returns Err rather
    // than leaving sibling threads stuck in retry loops on a never-draining
    // queue.
    let stage_factory = Arc::new(stage_factory);
    let mut worker_handles = Vec::with_capacity(cfg.worker_threads);
    for wid in 0..cfg.worker_threads {
        let q_in_worker = q_in.clone();
        let q_out_worker = q_out.clone();
        let cancel_worker = cancel.clone();
        let factory = stage_factory.clone();
        let handle = thread::Builder::new()
            .name(format!("pipeline_worker_{wid}"))
            .spawn(move || {
                let cancel_for_catch = cancel_worker.clone();
                run_and_catch_panic(&cancel_for_catch, move || {
                    let stage = factory();
                    run_worker(stage, q_in_worker, q_out_worker, cancel_worker)
                })
            })
            .context("failed to spawn worker thread")?;
        worker_handles.push(handle);
    }

    // Sink thread. Wraps the typed queue in a `CountingInputQueue` so every
    // successful pop increments the consumed counter.
    let cancel_for_sink = cancel.clone();
    let q_out_for_close = q_out.clone();
    let sink_in: Box<dyn InputQueue<Sk::Input>> = Box::new(CountingInputQueue::new(
        Box::new(q_out) as Box<dyn InputQueue<Sk::Input>>,
        consumed.clone(),
    ));
    let sink_handle = thread::Builder::new()
        .name("pipeline_sink".into())
        .spawn(move || sink.run(sink_in, cancel_for_sink))
        .context("failed to spawn sink thread")?;

    // Join all threads, preserving the FIRST error observed. On any error or
    // panic, cancel the pipeline so sibling threads exit their retry loops.
    let mut first_error: Option<anyhow::Error> = None;

    if let Err(e) = join_thread_or_cancel(source_handle, cancel) {
        first_error.get_or_insert(e);
    }

    for h in worker_handles {
        if let Err(e) = join_thread_or_cancel(h, cancel) {
            first_error.get_or_insert(e);
        }
    }

    // Close q_out now that all workers have finished. Workers no longer close
    // the output themselves; closing here avoids a race where the first worker
    // to finish would close q_out while a sibling still had an item in flight.
    q_out_for_close.close();

    if let Err(e) = join_thread_or_cancel(sink_handle, cancel) {
        first_error.get_or_insert(e);
    }

    // All queues have been manually closed above; nothing left for the guard.
    close_guard.disarm();

    if let Some(e) = first_error {
        return Err(e);
    }

    // Item-count equality (source produced == sink consumed) is NOT a valid
    // invariant once stages can coalesce or suppress outputs (e.g., FastqPair
    // pairs N per-stream chunks into M ≤ N/2 templates). Completeness is
    // validated by the sink's own reorder-buffer drain check instead.
    let _ = (&produced, &consumed, &cancel);

    // NOTE: This is the peak bytes held inside inter-stage queues, summed via
    // each item's `Stage::output_memory_estimate`. It is NOT the process RSS:
    // payloads held in sink / coalesce reorder buffers, worker held-item slots,
    // in-flight items currently being processed, stage-internal scratch
    // buffers, and subprocess memory (e.g. bwa mem) are not visible to the
    // tracker. Expect actual process RSS to be meaningfully higher.
    tracing::info!(
        "pipeline peak queued bytes: {} MiB (queue-held only; not process RSS)",
        global_memory.peak() / 1024 / 1024,
    );

    Ok(())
}

/// Type-erased stage interface. Allows heterogeneous stages to live in a single
/// `Vec<Box<dyn ErasedStage>>` so the driver can chain arbitrary numbers of
/// stages with differing `Input`/`Output` types.
pub trait ErasedStage: Send {
    /// Process one item. Input and output are type-erased through `Box<dyn Any>`.
    ///
    /// Calls `emit` with the output if the stage produces one. The callback
    /// may be called 0 or 1 times; calling it more than once panics in the
    /// [`TypedStage`] wrapper.
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying stage fails or if the input does not
    /// downcast to the expected type.
    fn process_any(
        &mut self,
        input: Box<dyn Any + Send>,
        emit: &mut dyn FnMut(Box<dyn Any + Send>),
    ) -> Result<()>;

    /// Memory estimate for a type-erased output item.
    fn memory_estimate_any(&self, output: &(dyn Any + Send)) -> usize;

    /// Parallelism constraint for this stage.
    fn parallelism(&self) -> super::stage::Parallelism;

    /// Human-readable stage name for logs and diagnostics.
    fn name(&self) -> &'static str;
}

/// Adapter that wraps any concrete [`Stage`] as an [`ErasedStage`].
///
/// Holds the inner stage plus `PhantomData` markers for its `Input`/`Output`
/// types so callers don't need to name them after construction.
pub struct TypedStage<S: Stage> {
    inner: S,
    _phantom: PhantomData<fn(S::Input) -> S::Output>,
}

impl<S: Stage> TypedStage<S> {
    /// Wrap a concrete stage.
    pub fn new(inner: S) -> Self {
        Self { inner, _phantom: PhantomData }
    }
}

impl<S: Stage + 'static> ErasedStage for TypedStage<S> {
    fn process_any(
        &mut self,
        input: Box<dyn Any + Send>,
        emit: &mut dyn FnMut(Box<dyn Any + Send>),
    ) -> Result<()> {
        let typed_input: Box<S::Input> = input.downcast::<S::Input>().map_err(|_| {
            anyhow::anyhow!("TypedStage '{}' received input of wrong type", self.inner.name())
        })?;
        self.inner.process(*typed_input, &mut |output| {
            emit(Box::new(output));
        })?;
        Ok(())
    }

    fn memory_estimate_any(&self, output: &(dyn Any + Send)) -> usize {
        match output.downcast_ref::<S::Output>() {
            Some(typed_output) => self.inner.output_memory_estimate(typed_output),
            None => std::mem::size_of::<S::Output>(),
        }
    }

    fn parallelism(&self) -> super::stage::Parallelism {
        self.inner.parallelism()
    }

    fn name(&self) -> &'static str {
        self.inner.name()
    }
}

/// Queue type carrying type-erased items between stages.
pub type ErasedQueue = StageQueue<Box<dyn Any + Send>>;

/// A stage entry in the pipeline chain.
///
/// Most stages are pool-managed ([`StageEntry::Pool`]): the work-stealing pool
/// assigns items to workers one at a time via [`ErasedStageSource`]. A
/// [`StageEntry::Special`] stage owns its execution loop and is given a
/// dedicated thread with direct queue handles.
pub enum StageEntry {
    /// A stage managed by the work-stealing pool.
    Pool(ErasedStageSource),
    /// A self-threaded stage with dedicated I/O queues.
    Special(Box<dyn ErasedSpecialStage>),
}

impl StageEntry {
    /// Human-readable name for logging.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            StageEntry::Pool(src) => src.name(),
            StageEntry::Special(s) => s.name(),
        }
    }
}

/// A contiguous segment of the pipeline chain.
enum Segment {
    /// Contiguous run of pool-managed stages.
    Pool(Vec<ErasedStageSource>),
    /// A single self-threaded special stage.
    Special(Box<dyn ErasedSpecialStage>),
}

/// Partition a flat list of [`StageEntry`] into contiguous segments.
///
/// Adjacent `Pool` entries are merged into a single `Segment::Pool`. Each
/// `Special` entry becomes its own `Segment::Special`. This lets the engine
/// create one work-stealing pool per contiguous pool-segment and spawn a
/// dedicated thread for each special stage.
fn partition_stages(entries: Vec<StageEntry>) -> Vec<Segment> {
    let mut segments: Vec<Segment> = Vec::new();
    let mut current_pool: Vec<ErasedStageSource> = Vec::new();

    for entry in entries {
        match entry {
            StageEntry::Pool(src) => current_pool.push(src),
            StageEntry::Special(special) => {
                if !current_pool.is_empty() {
                    segments.push(Segment::Pool(std::mem::take(&mut current_pool)));
                }
                segments.push(Segment::Special(special));
            }
        }
    }
    if !current_pool.is_empty() {
        segments.push(Segment::Pool(current_pool));
    }
    segments
}

/// Run a multi-stage pipeline: Source -> Stage1 -> ... -> `StageN` -> Sink
/// via the work-stealing pool.
///
/// Spawns:
/// - One source thread that pushes directly into the first erased inter-stage
///   queue via an [`ErasingOutputQueue`] wrapper. Items box exactly once here.
/// - A single persistent pool of `cfg.worker_threads` workers. Parallel
///   stages give each worker a per-worker instance (Clone or factory) via
///   [`ErasedStageSource`]; Sequential stages share a single Mutex-wrapped
///   instance so at most one worker enters them at a time.
/// - One sink thread that pops from the last erased inter-stage queue via an
///   [`UnerasingInputQueue`] wrapper. Items unbox exactly once here.
///
/// The driver thread owns closing every intermediate queue. Producer threads
/// (source, workers) never close their output queue themselves — this
/// prevents the race where a fast producer closes before a sibling finishes
/// pushing.
///
/// # Errors
///
/// Returns an error if any thread fails or panics, or if `stages` is empty.
pub fn run_multi_stage<SourceOut, SinkIn>(
    source: Box<dyn Source<Output = SourceOut>>,
    stages: Vec<StageEntry>,
    sink: Box<dyn Sink<Input = SinkIn>>,
    cancel: &CancelToken,
    cfg: PipelineConfig,
) -> Result<()>
where
    SourceOut: Send + 'static,
    SinkIn: Send + 'static,
{
    run_multi_stage_inner::<SourceOut, SinkIn>(source, stages, sink, cancel, cfg, None, None)
}

/// Variant of [`run_multi_stage`] that also writes per-stage metrics to a
/// TSV file after the pipeline completes successfully.
///
/// Columns: `stage`, `wall_time_secs`, `records_in`, `records_out`. One
/// row is emitted per pool stage in chain order; special (self-threaded)
/// stages are not instrumented and are therefore not emitted.
///
/// # Errors
///
/// Returns an error if the pipeline fails, a worker panics, or the TSV
/// cannot be written.
pub fn run_multi_stage_with_metrics<SourceOut, SinkIn>(
    source: Box<dyn Source<Output = SourceOut>>,
    stages: Vec<StageEntry>,
    sink: Box<dyn Sink<Input = SinkIn>>,
    cancel: &CancelToken,
    cfg: PipelineConfig,
    metrics_tsv: Option<&Path>,
) -> Result<()>
where
    SourceOut: Send + 'static,
    SinkIn: Send + 'static,
{
    run_multi_stage_inner::<SourceOut, SinkIn>(source, stages, sink, cancel, cfg, None, metrics_tsv)
}

/// Like [`run_multi_stage_with_metrics`] but additionally spawns a deadlock
/// watchdog. If no stage makes progress within `wd_cfg.timeout`, the pipeline
/// is cancelled and this function returns an error mentioning "deadlock".
///
/// # Errors
/// Returns an error if the pipeline fails, a worker panics, the TSV cannot
/// be written, or the watchdog trips.
pub fn run_multi_stage_with_metrics_and_watchdog<SourceOut, SinkIn>(
    source: Box<dyn Source<Output = SourceOut>>,
    stages: Vec<StageEntry>,
    sink: Box<dyn Sink<Input = SinkIn>>,
    cancel: &CancelToken,
    cfg: PipelineConfig,
    metrics_tsv: Option<&Path>,
    watchdog_cfg: WatchdogConfig,
) -> Result<()>
where
    SourceOut: Send + 'static,
    SinkIn: Send + 'static,
{
    let n_queues = stages.len() + 1;
    let progress: Vec<Arc<QueueProgress>> =
        (0..n_queues).map(|_| Arc::new(QueueProgress::default())).collect();

    let watchdog_exit = CancelToken::new();
    let watchdog_exit_for_thread = watchdog_exit.clone();
    let external_cancel_for_deadlock = cancel.clone();
    let watchdog_progress = progress.clone();
    let deadlock_flag = Arc::new(std::sync::atomic::AtomicBool::new(false));
    let deadlock_flag_for_watchdog = deadlock_flag.clone();

    let watchdog_handle = thread::Builder::new()
        .name("pipeline_watchdog".into())
        .spawn(move || {
            let outcome = watch(&watchdog_progress, &watchdog_exit_for_thread, watchdog_cfg);
            if outcome == WatchdogOutcome::Deadlock {
                deadlock_flag_for_watchdog.store(true, std::sync::atomic::Ordering::SeqCst);
                external_cancel_for_deadlock.cancel();
            }
            outcome
        })
        .context("failed to spawn watchdog thread")?;

    let run_result =
        run_multi_stage_inner(source, stages, sink, cancel, cfg, Some(&progress), metrics_tsv);

    watchdog_exit.cancel();
    let _ = watchdog_handle.join();

    if deadlock_flag.load(std::sync::atomic::Ordering::SeqCst) {
        anyhow::bail!(
            "pipeline deadlock detected: no progress for {}s",
            watchdog_cfg.timeout.as_secs()
        );
    }
    run_result
}

/// A running segment of the pipeline: either a pool (with its out-queue index)
/// or a special-stage thread (with its out-queue index). Used by
/// `run_multi_stage_inner` to join in chain order.
enum SpawnedSegment {
    Pool { pool: Pool, out_queue_idx: usize },
    Special { handle: thread::JoinHandle<Result<()>>, out_queue_idx: usize },
}

/// Internal implementation: wires the source and sink through the pool
/// driver, with optional per-queue progress counters used by the
/// deadlock-watchdog.
///
/// Expected `progress` length is `stages.len() + 1` — one counter per
/// erased inter-stage queue (input of stage 0, output of stage 0 / input
/// of stage 1, ..., output of the final stage / input of the sink).
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn run_multi_stage_inner<SourceOut, SinkIn>(
    source: Box<dyn Source<Output = SourceOut>>,
    stages: Vec<StageEntry>,
    sink: Box<dyn Sink<Input = SinkIn>>,
    cancel: &CancelToken,
    cfg: PipelineConfig,
    progress: Option<&Vec<Arc<QueueProgress>>>,
    metrics_tsv: Option<&Path>,
) -> Result<()>
where
    SourceOut: Send + 'static,
    SinkIn: Send + 'static,
{
    anyhow::ensure!(!stages.is_empty(), "run_multi_stage requires at least one stage");

    let global_memory = Arc::new(MemoryTracker::new(cfg.global_memory_limit));
    let n = stages.len();

    // Collect pool-stage names (special stages don't have pool workers and
    // aren't tracked by `PipelineStats`). This aligned list is used both for
    // the `format_summary` log emission and the `--metrics` TSV output.
    let pool_stage_names: Vec<&'static str> = stages
        .iter()
        .filter_map(|e| if matches!(e, StageEntry::Pool(_)) { Some(e.name()) } else { None })
        .collect();

    // Count only pool stages for PipelineStats (special stages have no pool workers).
    let num_pool_stages: usize = stages.iter().filter(|e| matches!(e, StageEntry::Pool(_))).count();

    // Create pipeline stats: num_queues = n + 1 (one queue per stage pair).
    let num_queues = n + 1;
    let mut stats_owned = PipelineStats::new(num_pool_stages, cfg.worker_threads, num_queues);
    // Register per-stage names so the pool can tag progress events by name.
    stats_owned.set_stage_names(pool_stage_names.iter().map(|s| (*s).to_string()).collect());
    let stats = Arc::new(stats_owned);

    // Announce each pool stage to the progress tracker. Ownership of
    // rendering (dashboard / heartbeat / no-op) is decided by the tracker
    // thread based on the `Mode` selected by `progress::init` upstream; this
    // driver only emits events.
    for name in &pool_stage_names {
        crate::progress::stage_started(name, None);
    }

    if let Some(p) = &progress {
        anyhow::ensure!(
            p.len() == n + 1,
            "progress vec length must be stages.len() + 1 (got {}, expected {})",
            p.len(),
            n + 1,
        );
    }

    // Erased inter-stage queues: `n + 1` total. queue[0] feeds stage 0 (from
    // the source); queue[i+1] carries stage i's output to stage i+1 (or the
    // sink for i == n-1).
    let mut erased_queues: Vec<Arc<ErasedQueue>> = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let name: String =
            if i == 0 { "stage_0_in".to_string() } else { format!("stage_{}_out", i - 1) };
        let q: Arc<ErasedQueue> = match &progress {
            Some(p) => Arc::new(StageQueue::with_progress(
                name,
                cfg.queue_capacity,
                cfg.queue_memory_limit,
                global_memory.clone(),
                p[i].clone(),
            )),
            None => Arc::new(StageQueue::new(
                name,
                cfg.queue_capacity,
                cfg.queue_memory_limit,
                global_memory.clone(),
            )),
        };
        erased_queues.push(q);
    }

    // Close-on-drop guard so a panic or early-return anywhere below still
    // unblocks every retry loop.
    let mut close_guard = QueueCloseGuard::new();
    for q in &erased_queues {
        close_guard.register(q.clone());
    }

    // Source/sink item count validation.
    let produced = Arc::new(AtomicU64::new(0));
    let consumed = Arc::new(AtomicU64::new(0));

    // === Source thread ===
    // Registers its well-known name so gate tests can enumerate pipeline
    // threads. Pushes into the first erased queue via an ErasingOutputQueue.
    let erasing_output: Box<dyn OutputQueue<SourceOut>> =
        Box::new(ErasingOutputQueue::<SourceOut>::new(erased_queues[0].clone()));
    let source_out: Box<dyn OutputQueue<SourceOut>> =
        Box::new(CountingOutputQueue::new(erasing_output, produced.clone()));
    let cancel_for_source = cancel.clone();
    let source_handle = thread::Builder::new()
        .name("pipeline_source".into())
        .spawn(move || -> Result<()> {
            test_support::register_thread("pipeline_source");
            source.run(source_out, cancel_for_source)
        })
        .context("failed to spawn source thread")?;

    // === Partition stages into segments ===
    // Each segment is either a contiguous run of pool stages or a single
    // special stage. Pool segments get their own work-stealing pool; special
    // stages run on dedicated threads.
    let segments = partition_stages(stages);

    // Count pool segments for proportional thread-budget splitting.
    let total_pool_stages: usize = segments
        .iter()
        .map(|seg| match seg {
            Segment::Pool(v) => v.len(),
            Segment::Special(_) => 0,
        })
        .sum();

    // Assign thread budgets to pool segments proportionally by stage count,
    // minimum 1 each.
    let num_pool_segs: usize = segments.iter().filter(|s| matches!(s, Segment::Pool(_))).count();
    let mut pool_thread_budgets: Vec<usize> = Vec::new();
    if num_pool_segs > 0 {
        let mut remaining_threads = cfg.worker_threads;
        let pool_stage_counts: Vec<usize> = segments
            .iter()
            .filter_map(|s| if let Segment::Pool(v) = s { Some(v.len()) } else { None })
            .collect();
        for (idx, &stage_count) in pool_stage_counts.iter().enumerate() {
            let threads = if idx == num_pool_segs - 1 {
                // Last pool segment gets whatever threads remain (minimum 1).
                remaining_threads.max(1)
            } else {
                // Proportional share, clamped so at least 1 thread remains for
                // the segments that follow.
                let share = (cfg.worker_threads * stage_count / total_pool_stages).max(1);
                share.min(remaining_threads.saturating_sub(1).max(1))
            };
            pool_thread_budgets.push(threads);
            remaining_threads = remaining_threads.saturating_sub(threads).max(1);
        }
    }

    // === Spawn segment threads ===
    // We track: a list of (pool, queue_range) and special thread handles.
    // After source joins, we close queue[0], then join segments in order,
    // closing the output queue of each segment before moving to the next.
    let mut spawned: Vec<SpawnedSegment> = Vec::new();
    let mut pool_budget_iter = pool_thread_budgets.into_iter();

    // Track the queue index consumed by each segment's input side.
    // Segment i reads from erased_queues[i] and writes to erased_queues[i+1].
    let mut seg_queue_start: usize = 0; // index of the first queue for next segment

    for segment in segments {
        match segment {
            Segment::Pool(pool_stages) => {
                let seg_n = pool_stages.len();
                let thread_budget = pool_budget_iter.next().unwrap_or(cfg.worker_threads);

                // The queues for this pool segment: [seg_queue_start .. seg_queue_start + seg_n + 1].
                // queue[seg_queue_start] is the input; queue[seg_queue_start + seg_n] is the output.
                let seg_queues = erased_queues[seg_queue_start..=seg_queue_start + seg_n].to_vec();
                let out_queue_idx = seg_queue_start + seg_n;

                let invokers: Vec<StageInvoker> =
                    pool_stages.into_iter().map(StageInvoker::from_source).collect();
                let pool = Pool::start(
                    invokers,
                    seg_queues,
                    cancel,
                    thread_budget,
                    global_memory.clone(),
                    &stats,
                )?;

                spawned.push(SpawnedSegment::Pool { pool, out_queue_idx });
                seg_queue_start += seg_n;
            }
            Segment::Special(special) => {
                let in_q = erased_queues[seg_queue_start].clone();
                let out_q = erased_queues[seg_queue_start + 1].clone();
                let out_queue_idx = seg_queue_start + 1;
                let cancel_for_special = cancel.clone();
                let stage_name: &'static str = special.name();

                let handle = thread::Builder::new()
                    .name(format!("pipeline_special_{stage_name}"))
                    .spawn(move || -> Result<()> {
                        run_and_catch_panic_named(
                            &cancel_for_special.clone(),
                            Some(stage_name),
                            move || special.run_erased(in_q, out_q, cancel_for_special),
                        )
                    })
                    .context("failed to spawn special stage thread")?;

                spawned.push(SpawnedSegment::Special { handle, out_queue_idx });
                seg_queue_start += 1;
            }
        }
    }

    // === Sink thread ===
    let unerasing_input: Box<dyn InputQueue<SinkIn>> =
        Box::new(UnerasingInputQueue::<SinkIn>::new(erased_queues[n].clone()));
    let sink_in: Box<dyn InputQueue<SinkIn>> =
        Box::new(CountingInputQueue::new(unerasing_input, consumed.clone()));
    let cancel_for_sink = cancel.clone();
    let sink_handle = thread::Builder::new()
        .name("pipeline_sink".into())
        .spawn(move || -> Result<()> {
            test_support::register_thread("pipeline_sink");
            sink.run(sink_in, cancel_for_sink)
        })
        .context("failed to spawn sink thread")?;

    // === Join: source → segments in order → sink. ===
    // Propagate the first error; cancel the pipeline on any failure so
    // siblings exit their retry loops.
    let mut first_error: Option<anyhow::Error> = None;

    // 1. Source: closes queue[0] on its own per the Source contract.
    if let Err(e) = join_thread_or_cancel(source_handle, cancel) {
        first_error.get_or_insert(e);
    }
    erased_queues[0].close();

    // 2. Join segments in chain order. After each segment finishes, close its
    //    output queue so the next segment (or sink) notices drain.
    for seg in spawned {
        match seg {
            SpawnedSegment::Pool { pool, out_queue_idx } => {
                if let Err(e) = pool.join_all(cancel) {
                    first_error.get_or_insert(e);
                }
                erased_queues[out_queue_idx].close();
            }
            SpawnedSegment::Special { handle, out_queue_idx } => {
                // Special stages call output.close() themselves per the
                // SpecialStage contract, but we close here too to be safe
                // in error paths (idempotent).
                if let Err(e) = join_thread_or_cancel(handle, cancel) {
                    first_error.get_or_insert(e);
                }
                erased_queues[out_queue_idx].close();
            }
        }
    }

    // 3. Sink.
    if let Err(e) = join_thread_or_cancel(sink_handle, cancel) {
        first_error.get_or_insert(e);
    }

    close_guard.disarm();

    // Emit stage_finished events for every pool stage so the tracker can
    // close its bars (dashboard) or summary rows (heartbeat). Done regardless
    // of error so a cancelled/failed pipeline still lands a clean tracker
    // state for the summary.
    for name in &pool_stage_names {
        crate::progress::stage_finished(name);
    }

    if let Some(e) = first_error {
        return Err(e);
    }

    // Item-count equality is NOT a valid invariant once stages can coalesce
    // or suppress outputs; see the non-watchdog variant above. Completeness
    // is validated by the sink's own reorder-buffer drain check.
    let _ = (&produced, &consumed, &cancel);

    tracing::info!("{}", stats.format_summary(&pool_stage_names));
    // See accompanying note in `run_multi_stage` — this tracks bytes held in
    // queues, not process RSS.
    tracing::info!(
        "pipeline peak queued bytes: {} MiB (queue-held only; not process RSS)",
        global_memory.peak() / 1024 / 1024,
    );

    // Write per-stage metrics TSV if requested.
    if let Some(path) = metrics_tsv {
        stats
            .write_tsv(path, &pool_stage_names)
            .with_context(|| format!("failed to write metrics TSV to {}", path.display()))?;
    }

    Ok(())
}

/// Variant of [`run_multi_stage`] with optional watchdog.
///
/// When `watchdog_cfg` is `Some`, spawns a watchdog thread that monitors
/// progress. If the watchdog times out with no forward progress, the pipeline
/// is cancelled and this function returns an error mentioning "deadlock".
///
/// # Errors
/// Returns an error if the pipeline fails or a deadlock is detected.
pub fn run_multi_stage_with_watchdog<SourceOut, SinkIn>(
    source: Box<dyn Source<Output = SourceOut>>,
    stages: Vec<StageEntry>,
    sink: Box<dyn Sink<Input = SinkIn>>,
    cancel: &CancelToken,
    cfg: PipelineConfig,
    watchdog_cfg: Option<WatchdogConfig>,
) -> Result<()>
where
    SourceOut: Send + 'static,
    SinkIn: Send + 'static,
{
    let Some(wd_cfg) = watchdog_cfg else {
        return run_multi_stage(source, stages, sink, cancel, cfg);
    };

    // One QueueProgress per erased inter-stage queue. The
    // run_multi_stage_inner entry point attaches each counter to its
    // matching queue, so the watchdog sees real push/pop traffic instead of
    // a static placeholder. We allocate the counters here and pass them in.
    let n_queues = stages.len() + 1; // n + 1 erased queues (source-in .. sink-in)
    let progress: Vec<Arc<QueueProgress>> =
        (0..n_queues).map(|_| Arc::new(QueueProgress::default())).collect();

    // Two separate tokens govern the watchdog's lifecycle:
    //
    // - `watchdog_exit` is internal to this function. We cancel it when the
    //   pipeline has returned, to signal the watchdog thread to exit. It must
    //   NOT affect the caller's cancel token.
    // - `external_cancel_for_deadlock` is a clone of the caller's token. The
    //   watchdog cancels this ONLY when it detects a deadlock — so callers
    //   never observe a spurious cancellation on a successful pipeline run.
    let watchdog_exit = CancelToken::new();
    let watchdog_exit_for_thread = watchdog_exit.clone();
    let external_cancel_for_deadlock = cancel.clone();
    let watchdog_progress = progress.clone();
    let deadlock_flag = Arc::new(std::sync::atomic::AtomicBool::new(false));
    let deadlock_flag_for_watchdog = deadlock_flag.clone();

    let watchdog_handle = thread::Builder::new()
        .name("pipeline_watchdog".into())
        .spawn(move || {
            let outcome = watch(&watchdog_progress, &watchdog_exit_for_thread, wd_cfg);
            if outcome == WatchdogOutcome::Deadlock {
                deadlock_flag_for_watchdog.store(true, std::sync::atomic::Ordering::SeqCst);
                external_cancel_for_deadlock.cancel();
            }
            outcome
        })
        .context("failed to spawn watchdog thread")?;

    let run_result =
        run_multi_stage_inner(source, stages, sink, cancel, cfg, Some(&progress), None);

    // Tell the watchdog to exit. This touches only our internal token — the
    // caller's `cancel` is unaffected on a successful run.
    watchdog_exit.cancel();
    let _ = watchdog_handle.join();

    if deadlock_flag.load(std::sync::atomic::Ordering::SeqCst) {
        anyhow::bail!("pipeline deadlock detected: no progress for {}s", wd_cfg.timeout.as_secs());
    }
    run_result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::stage::{Parallelism, SequencedItem};
    use std::sync::Mutex;

    #[test]
    fn test_pipeline_config_defaults() {
        let cfg = PipelineConfig::default();
        assert_eq!(cfg.worker_threads, 1);
        assert_eq!(cfg.queue_capacity, DEFAULT_QUEUE_CAPACITY);
        assert_eq!(cfg.queue_capacity, 1024);
        assert_eq!(cfg.queue_memory_limit, DEFAULT_QUEUE_MEMORY_LIMIT);
        assert_eq!(cfg.queue_memory_limit, 2 * 1024 * 1024 * 1024);
        assert_eq!(cfg.global_memory_limit, DEFAULT_GLOBAL_MEMORY_LIMIT);
        assert_eq!(cfg.global_memory_limit, 4 * 1024 * 1024 * 1024);
    }

    #[test]
    fn test_auto_tuned_queue_capacity_clamps_at_low_end() {
        // 1 thread -> 1*16 = 16, clamped up to 64.
        assert_eq!(PipelineConfig::auto_tuned(1).queue_capacity, 64);
        assert_eq!(PipelineConfig::auto_tuned(1).worker_threads, 1);
        assert_eq!(PipelineConfig::auto_tuned(4).queue_capacity, 64);
        // 2 threads -> 2*16 = 32, still clamped up to 64.
        assert_eq!(PipelineConfig::auto_tuned(2).queue_capacity, 64);
    }

    #[test]
    fn test_auto_tuned_queue_capacity_scales_in_the_middle() {
        assert_eq!(PipelineConfig::auto_tuned(8).queue_capacity, 128);
        assert_eq!(PipelineConfig::auto_tuned(10).queue_capacity, 160);
    }

    #[test]
    fn test_auto_tuned_queue_capacity_saturates_at_high_end() {
        assert_eq!(PipelineConfig::auto_tuned(16).queue_capacity, 256);
        assert_eq!(PipelineConfig::auto_tuned(32).queue_capacity, 256);
        assert_eq!(PipelineConfig::auto_tuned(1024).queue_capacity, 256);
    }

    #[test]
    fn test_auto_tuned_preserves_default_memory_limits() {
        let cfg = PipelineConfig::auto_tuned(8);
        assert_eq!(cfg.queue_memory_limit, DEFAULT_QUEUE_MEMORY_LIMIT);
        assert_eq!(cfg.global_memory_limit, DEFAULT_GLOBAL_MEMORY_LIMIT);
        assert_eq!(cfg.worker_threads, 8);
    }

    struct RangeSource {
        n: u64,
    }

    impl Source for RangeSource {
        type Output = u64;

        fn run(
            self: Box<Self>,
            output: Box<dyn super::super::source::OutputQueue<Self::Output>>,
            cancel: CancelToken,
        ) -> anyhow::Result<()> {
            for seq in 0..self.n {
                if cancel.is_cancelled() {
                    break;
                }
                let item = SequencedItem::new(seq, seq, 8);
                if output.push_until_cancelled(item, &cancel).is_err() {
                    break;
                }
            }
            output.close();
            Ok(())
        }

        fn name(&self) -> &'static str {
            "RangeSource"
        }
    }

    #[derive(Clone)]
    struct SquareStage;

    impl Stage for SquareStage {
        type Input = u64;
        type Output = u64;

        fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
            out(input * input);
            Ok(())
        }

        fn parallelism(&self) -> Parallelism {
            Parallelism::Parallel
        }

        fn output_memory_estimate(&self, _output: &u64) -> usize {
            std::mem::size_of::<u64>()
        }

        fn name(&self) -> &'static str {
            "SquareStage"
        }
    }

    struct CollectSink {
        out: Arc<Mutex<Vec<u64>>>,
    }

    impl Sink for CollectSink {
        type Input = u64;

        fn run(
            self: Box<Self>,
            input: Box<dyn super::super::sink::InputQueue<Self::Input>>,
            cancel: CancelToken,
        ) -> anyhow::Result<()> {
            let mut buffer: Vec<(u64, u64)> = Vec::new();
            loop {
                if cancel.is_cancelled() {
                    break;
                }
                if let Some(item) = input.pop() {
                    buffer.push((item.seq, item.item));
                } else {
                    if input.is_drained() {
                        break;
                    }
                    std::thread::yield_now();
                }
            }
            buffer.sort_by_key(|(seq, _)| *seq);
            let mut out = self.out.lock().unwrap();
            for (_, v) in buffer {
                out.push(v);
            }
            Ok(())
        }

        fn name(&self) -> &'static str {
            "CollectSink"
        }
    }

    #[test]
    fn test_single_stage_end_to_end() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out: out.clone() });
        let source = Box::new(RangeSource { n: 10 });

        run_single_stage(
            source,
            || SquareStage,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 2,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 1_000_000,
            },
        )
        .unwrap();

        let got = out.lock().unwrap().clone();
        let expected: Vec<_> = (0..10u64).map(|i| i * i).collect();
        assert_eq!(got, expected);
    }

    #[derive(Clone)]
    struct IncrementStage;

    impl Stage for IncrementStage {
        type Input = u64;
        type Output = u64;

        fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
            out(input + 1);
            Ok(())
        }

        fn parallelism(&self) -> Parallelism {
            Parallelism::Parallel
        }

        fn output_memory_estimate(&self, _output: &u64) -> usize {
            std::mem::size_of::<u64>()
        }

        fn name(&self) -> &'static str {
            "IncrementStage"
        }
    }

    #[test]
    fn test_two_stage_pipeline() {
        // Source: 0..10. Stage 1: square. Stage 2: +1. Sink: collect.
        // Expected: [1, 2, 5, 10, 17, 26, 37, 50, 65, 82]
        use crate::runall::engine::stage_source::{ErasedStageSource, StageSource};
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out: out.clone() });
        let source = Box::new(RangeSource { n: 10 });

        let stages: Vec<StageEntry> = vec![
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(SquareStage))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(IncrementStage))),
        ];

        run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 2,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 1_000_000,
            },
        )
        .unwrap();

        let got = out.lock().unwrap().clone();
        let expected: Vec<_> = (0..10u64).map(|i| i * i + 1).collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn test_cancellation_stops_pipeline() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out: out.clone() });
        let source = Box::new(RangeSource { n: 1_000_000 });
        let cancel = CancelToken::new();

        // Cancel shortly after starting.
        let cancel_for_timer = cancel.clone();
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(10));
            cancel_for_timer.cancel();
        });

        run_single_stage(source, || SquareStage, sink, &cancel, PipelineConfig::default()).unwrap();

        // We don't assert on count — just that we finish promptly (not 1M items).
        let count = out.lock().unwrap().len();
        assert!(count < 1_000_000);
    }

    #[derive(Clone)]
    struct BlockingStage {
        block_until_cancel: std::sync::Arc<std::sync::atomic::AtomicBool>,
    }

    impl Stage for BlockingStage {
        type Input = u64;
        type Output = u64;

        fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
            while !self.block_until_cancel.load(std::sync::atomic::Ordering::Relaxed) {
                std::thread::sleep(std::time::Duration::from_millis(10));
            }
            out(input);
            Ok(())
        }

        fn parallelism(&self) -> Parallelism {
            Parallelism::Parallel
        }

        fn output_memory_estimate(&self, _output: &u64) -> usize {
            std::mem::size_of::<u64>()
        }

        fn name(&self) -> &'static str {
            "BlockingStage"
        }
    }

    #[test]
    fn test_watchdog_detects_deadlocked_pipeline() {
        use crate::runall::engine::stage_source::{ErasedStageSource, StageSource};
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out: out.clone() });
        let source = Box::new(RangeSource { n: 10 });
        let block = Arc::new(std::sync::atomic::AtomicBool::new(false));

        let stages: Vec<StageEntry> = vec![StageEntry::Pool(ErasedStageSource::from_source(
            StageSource::Clone(BlockingStage { block_until_cancel: block.clone() }),
        ))];

        let cfg = PipelineConfig {
            worker_threads: 1,
            queue_capacity: 2,
            queue_memory_limit: 1024,
            global_memory_limit: 1_000_000,
        };

        let cancel = CancelToken::new();

        // Short watchdog timeout so the test finishes quickly.
        let watchdog_cfg = crate::runall::engine::deadlock::WatchdogConfig {
            timeout: std::time::Duration::from_millis(300),
            poll_interval: std::time::Duration::from_millis(30),
        };

        // BlockingStage::process() cannot observe the pipeline cancel — it loops
        // on its own flag. Start a background timer that unblocks the stage
        // AFTER the watchdog has had time to fire (300ms timeout + margin).
        // This lets the pipeline unwind cleanly once the deadlock is detected.
        let block_unblocker = block.clone();
        let unblocker = std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(600));
            block_unblocker.store(true, std::sync::atomic::Ordering::Relaxed);
        });

        let result = run_multi_stage_with_watchdog::<u64, u64>(
            source,
            stages,
            sink,
            &cancel,
            cfg,
            Some(watchdog_cfg),
        );

        unblocker.join().expect("unblocker thread panicked");

        // Expect a deadlock error, not a successful completion, and not a hang.
        assert!(result.is_err(), "expected deadlock error, got {result:?}");
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.to_lowercase().contains("deadlock"), "expected deadlock message, got: {msg}");
    }

    #[test]
    fn test_watchdog_sees_progress_from_queues() {
        // Drive a well-behaved pipeline with a short watchdog timeout but a
        // long enough poll interval that traffic will be observed between
        // polls. With real per-queue progress instrumentation, the watchdog
        // must NOT report a false-positive deadlock.
        use crate::runall::engine::stage_source::{ErasedStageSource, StageSource};
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out: out.clone() });
        let source = Box::new(RangeSource { n: 200 });

        let stages: Vec<StageEntry> = vec![StageEntry::Pool(ErasedStageSource::from_source(
            StageSource::Clone(IncrementStage),
        ))];

        let cfg = PipelineConfig {
            worker_threads: 1,
            queue_capacity: 4,
            queue_memory_limit: 1024,
            global_memory_limit: 1_000_000,
        };

        let cancel = CancelToken::new();
        let watchdog_cfg = crate::runall::engine::deadlock::WatchdogConfig {
            timeout: std::time::Duration::from_millis(500),
            poll_interval: std::time::Duration::from_millis(20),
        };

        let result = run_multi_stage_with_watchdog::<u64, u64>(
            source,
            stages,
            sink,
            &cancel,
            cfg,
            Some(watchdog_cfg),
        );

        assert!(result.is_ok(), "expected success, got {result:?}");
        assert_eq!(out.lock().unwrap().len(), 200);
    }
}

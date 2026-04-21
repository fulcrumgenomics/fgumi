//! Work-stealing thread pool for `pipeline`.
//!
//! Spawns `worker_threads` persistent OS threads. Each worker runs a loop:
//!
//! 1. Check cancellation.
//! 2. Sample backpressure across all inter-stage queues.
//! 3. Ask [`Scheduler::priorities_filtered`] for the order in which to try stages.
//! 4. For each stage in priority order, attempt:
//!    - flush any held item into the stage's output queue;
//!    - pop a new input, run `Stage::process`, push the output.
//! 5. On success, update affinity and reset backoff.
//! 6. On no progress anywhere, apply exponential backoff.
//!
//! Parallel stages give each worker a per-worker instance (via
//! [`ErasedStageSource::build_erased`]). Sequential stages are wrapped in a
//! `Mutex<Box<dyn ErasedStage>>`; workers use `try_lock` to skip the stage
//! when another worker holds the lock. Barrier and Special stages are NOT
//! run by this pool — the driver handles them on dedicated threads.

use std::collections::{HashMap, VecDeque};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use anyhow::Result;

use super::backoff::Backoff;
use super::backpressure::{BackpressureState, QueueSummary};
use super::cancel::CancelToken;
use super::driver::{ErasedQueue, ErasedStage, extract_panic_message, run_and_catch_panic};
use super::memory::MemoryTracker;
use super::scheduler::{Scheduler, StageAttemptOutcome, WorkerRole, assign_sequential_ownership};
use super::stage::{Parallelism, SequencedItem};
use super::stage_source::ErasedStageSource;
use super::stats::PipelineStats;
use super::test_support;

/// Snapshot each stage's `Parallelism` in order.
fn parallelism_snapshot(invokers: &[StageInvoker]) -> Vec<Parallelism> {
    invokers
        .iter()
        .map(|inv| match inv {
            StageInvoker::Parallel(_) => Parallelism::Parallel,
            StageInvoker::Sequential(_) => Parallelism::Sequential,
        })
        .collect()
}

/// Per-stage shared state: either a per-worker factory (Parallel) or a
/// `Mutex`-guarded single instance (Sequential).
pub enum StageInvoker {
    /// A Parallel stage: each worker holds its own instance constructed via
    /// the factory. The factory also produces new instances for workers that
    /// start mid-run (not currently used, but keeps the API symmetric).
    Parallel(Arc<ErasedStageSource>),
    /// A Sequential stage: a single mutex-guarded instance shared across the
    /// pool. At most one worker enters `process()` at a time.
    Sequential(Arc<Mutex<Box<dyn ErasedStage>>>),
}

impl StageInvoker {
    /// Build a [`StageInvoker`] from an [`ErasedStageSource`] according to the
    /// declared parallelism.
    ///
    /// # Panics
    /// Panics if `src.parallelism()` is `Barrier` or `Special` — those stage
    /// kinds live outside the pool and the driver must handle them on
    /// dedicated threads.
    #[must_use]
    pub fn from_source(src: ErasedStageSource) -> Self {
        match src.parallelism() {
            Parallelism::Parallel => StageInvoker::Parallel(Arc::new(src)),
            Parallelism::Sequential => {
                let instance = src.build_erased();
                StageInvoker::Sequential(Arc::new(Mutex::new(instance)))
            }
            Parallelism::Barrier | Parallelism::Special => {
                panic!(
                    "Barrier and Special stages must not be passed to the pool \
                     (stage '{}')",
                    src.name()
                );
            }
        }
    }
}

/// Per-worker state: its own scheduler, its own per-Parallel-stage instances,
/// its own held items, and its own backoff.
///
/// Note: `stage_exited` does NOT live here — it's owned by the spawning
/// thread outside the worker loop so that a panic unwinding through the
/// worker body doesn't invalidate it before the cleanup guard fires.
struct WorkerState {
    scheduler: Scheduler,
    /// Per-stage Parallel instance. `None` for Sequential stages — they live
    /// inside the shared `Mutex` in [`StageInvoker::Sequential`].
    parallel_instances: Vec<Option<Box<dyn ErasedStage>>>,
    /// Per-stage OUTPUT-held queue: items already PROCESSED but not yet pushed
    /// because the output queue was full. A stage may emit more than one item
    /// per `process()` call (e.g. `FastqPair` chunked emit); each item is
    /// queued here and drained front-to-back on subsequent iterations.
    held_output: HashMap<usize, VecDeque<SequencedItem<Box<dyn std::any::Any + Send>>>>,
    /// Per-stage INPUT-held slot: an item already POPPED from the input
    /// queue but not yet processed because a Sequential stage was
    /// contended. The worker retries on the next iteration.
    held_input: HashMap<usize, SequencedItem<Box<dyn std::any::Any + Send>>>,
    backoff: Backoff,
}

impl WorkerState {
    fn new(scheduler: Scheduler, invokers: &[StageInvoker]) -> Self {
        let parallel_instances: Vec<_> = invokers
            .iter()
            .map(|inv| match inv {
                StageInvoker::Parallel(src) => Some(src.build_erased()),
                StageInvoker::Sequential(_) => None,
            })
            .collect();
        Self {
            scheduler,
            parallel_instances,
            held_output: HashMap::new(),
            held_input: HashMap::new(),
            backoff: Backoff::new(),
        }
    }

    /// Attempt one iteration of stage `idx`. Returns whether any progress
    /// was made, or the kind of obstruction the worker ran into.
    ///
    /// Held-item ordering:
    /// 1. Flush the front of `held_output[idx]` first — completed items pending push.
    /// 2. Otherwise retry `held_input[idx]` — a popped-but-unprocessed item
    ///    held because a Sequential stage was contended.
    /// 3. Otherwise pop a fresh input and process it.
    fn try_run_stage(
        &mut self,
        idx: usize,
        invoker: &StageInvoker,
        input: &Arc<ErasedQueue>,
        output: &Arc<ErasedQueue>,
        cancel: &CancelToken,
        pool_stats: &PipelineStats,
    ) -> StageAttemptResult {
        if cancel.is_cancelled() {
            return StageAttemptResult::Cancelled;
        }

        // 1. Flush a pending output first. If the output queue is still full,
        //    put the item back and report backpressure. The held-output queue
        //    may have multiple items (from a multi-emit stage like FastqPair);
        //    drain one item per iteration to give the scheduler a chance to
        //    rebalance between iterations.
        if let Some(held_item) = self.held_output.get_mut(&idx).and_then(VecDeque::pop_front) {
            match output.push(held_item) {
                Ok(()) => {
                    if self.held_output.get(&idx).is_some_and(VecDeque::is_empty) {
                        self.held_output.remove(&idx);
                    }
                    pool_stats.record_stage_records_out(idx);
                    crate::progress::records_out(pool_stats.stage_name(idx), 1);
                    return StageAttemptResult::Progress;
                }
                Err(returned) => {
                    self.held_output.entry(idx).or_default().push_front(returned);
                    return StageAttemptResult::BackpressureOutput;
                }
            }
        }

        // 2. Retrieve a held input, or pop a fresh one.
        //
        //    We only count `records_in` on a FRESH pop (not on held-input
        //    retries) so the total matches the number of distinct items the
        //    stage has consumed from its input queue.
        let input_item = if let Some(item) = self.held_input.remove(&idx) {
            item
        } else {
            let Some(popped) = input.pop() else {
                if input.is_drained() {
                    return StageAttemptResult::Drained;
                }
                return StageAttemptResult::BackpressureInput;
            };
            pool_stats.record_stage_records_in(idx);
            popped
        };
        let seq = input_item.seq;

        // 3. Run the stage. Parallel stages use the per-worker instance;
        //    Sequential stages use try_lock on the shared mutex.
        //
        //    A stage may emit more than one output per process() call (e.g.
        //    FastqPair chunked emit). Collect all outputs; push the first
        //    immediately and stash any overflow in held_output for later
        //    iterations.
        let mut captured: Vec<Box<dyn std::any::Any + Send>> = Vec::new();

        match invoker {
            StageInvoker::Sequential(mutex) => {
                let Ok(mut guard) = mutex.try_lock() else {
                    // Contention: another worker holds the Sequential
                    // stage. Put the input into the input-held slot so we
                    // retry next iteration without re-popping.
                    self.held_input.insert(idx, input_item);
                    return StageAttemptResult::SequentialBusy;
                };
                match invoke_stage_catching_panic(&mut **guard, input_item.item, &mut captured) {
                    Ok(()) => {}
                    Err(e) => return StageAttemptResult::Error(e),
                }
                if captured.is_empty() {
                    // Suppress-on-empty: stage emitted nothing.
                    return StageAttemptResult::Progress;
                }
                // Build SequencedItems for all captured outputs, reusing seq
                // for all (downstream extract pipeline orders by TemplateBatch
                // ordinal, not by pool seq).
                let items: VecDeque<_> = captured
                    .into_iter()
                    .map(|v| {
                        let mem = guard.memory_estimate_any(v.as_ref());
                        SequencedItem::new(seq, v, mem)
                    })
                    .collect();
                drop(guard);
                // 4. Push the first output. On full output queue, stash all
                //    remaining (including the failed first) into held_output.
                self.push_or_hold(idx, items, output, pool_stats)
            }
            StageInvoker::Parallel(_) => {
                let stage = self.parallel_instances[idx]
                    .as_mut()
                    .expect("Parallel stage must have a per-worker instance");
                match invoke_stage_catching_panic(&mut **stage, input_item.item, &mut captured) {
                    Ok(()) => {}
                    Err(e) => return StageAttemptResult::Error(e),
                }
                if captured.is_empty() {
                    return StageAttemptResult::Progress;
                }
                let stage = self.parallel_instances[idx]
                    .as_ref()
                    .expect("Parallel stage must have a per-worker instance");
                let items: VecDeque<_> = captured
                    .into_iter()
                    .map(|v| {
                        let mem = stage.memory_estimate_any(v.as_ref());
                        SequencedItem::new(seq, v, mem)
                    })
                    .collect();
                self.push_or_hold(idx, items, output, pool_stats)
            }
        }
    }

    /// Push the first item in `items` to `output`; stash any overflow (plus
    /// the failed item on backpressure) into `held_output[idx]` for later
    /// iterations. Returns `Progress` or `BackpressureOutput`.
    fn push_or_hold(
        &mut self,
        idx: usize,
        mut items: VecDeque<SequencedItem<Box<dyn std::any::Any + Send>>>,
        output: &Arc<ErasedQueue>,
        pool_stats: &PipelineStats,
    ) -> StageAttemptResult {
        let first = items.pop_front().expect("items must be non-empty");
        match output.push(first) {
            Ok(()) => {
                pool_stats.record_stage_records_out(idx);
                crate::progress::records_out(pool_stats.stage_name(idx), 1);
                if !items.is_empty() {
                    self.held_output.insert(idx, items);
                }
                StageAttemptResult::Progress
            }
            Err(returned) => {
                items.push_front(returned);
                self.held_output.insert(idx, items);
                StageAttemptResult::BackpressureOutput
            }
        }
    }

    /// Whether any held items (input or output) remain in this worker's
    /// slots, meaning the worker still has unfinished work.
    fn has_held_items(&self) -> bool {
        !self.held_output.is_empty() || !self.held_input.is_empty()
    }

    /// Whether this worker holds any item for `stage_idx` that still needs
    /// pushing or processing.
    fn holds_for_stage(&self, stage_idx: usize) -> bool {
        self.held_output.contains_key(&stage_idx) || self.held_input.contains_key(&stage_idx)
    }
}

/// Invoke `stage.process_any` while catching any panic it raises, so a
/// broken stage body yields a stage-attributed `anyhow::Error` (routed
/// through the normal stage-error path) instead of unwinding out of the
/// pool worker and leaving siblings stuck.
///
/// On panic, the returned error is tagged with the stage name via
/// [`ErasedStage::name`]. On `Err(_)` from `process_any`, the error is
/// returned unchanged. On panic, [`std::panic::AssertUnwindSafe`] is
/// sound: the cancel token is set by the caller immediately afterwards,
/// and no shared state is re-accessed by siblings except through the
/// inherently `Sync` `StageQueue`.
fn invoke_stage_catching_panic(
    stage: &mut dyn ErasedStage,
    input: Box<dyn std::any::Any + Send>,
    captured: &mut Vec<Box<dyn std::any::Any + Send>>,
) -> Result<()> {
    let stage_name = stage.name();
    let process_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        stage.process_any(input, &mut |v| {
            captured.push(v);
        })
    }));
    match process_result {
        Ok(Ok(())) => Ok(()),
        Ok(Err(e)) => Err(e),
        Err(payload) => {
            let msg = extract_panic_message(payload.as_ref());
            Err(anyhow::anyhow!("stage '{stage_name}' panicked: {msg}"))
        }
    }
}

/// Result of one stage-attempt on a worker iteration.
enum StageAttemptResult {
    /// Made forward progress (pushed an output).
    Progress,
    /// Input queue was empty but not closed; retry later.
    BackpressureInput,
    /// Output queue was full; held the item for later.
    BackpressureOutput,
    /// Sequential stage was locked by another worker.
    SequentialBusy,
    /// Input queue is closed and empty; stage can do no more work.
    Drained,
    /// The cancel token is set.
    Cancelled,
    /// Stage returned an error.
    Error(anyhow::Error),
}

/// Drop guard that runs the cleanup-on-exit cascade for a single pool
/// worker: on every exit path (normal return, error return, panic unwind),
/// decrement `producers_remaining[i]` for every stage the worker has not
/// already exited, and close `queues[i+1]` if this worker was the last
/// producer for stage i. Guarantees downstream queues unblock even if the
/// worker panics mid-iteration.
struct ExitGuard {
    stage_exited: Arc<Mutex<Vec<bool>>>,
    producers_remaining: Arc<Vec<AtomicUsize>>,
    queues: Arc<Vec<Arc<ErasedQueue>>>,
}

impl Drop for ExitGuard {
    fn drop(&mut self) {
        let stage_exited = self.stage_exited.lock().unwrap();
        for (stage_idx, exited) in stage_exited.iter().enumerate() {
            if !*exited {
                let prev = self.producers_remaining[stage_idx].fetch_sub(1, Ordering::SeqCst);
                if prev == 1 {
                    self.queues[stage_idx + 1].close();
                }
            }
        }
    }
}

/// Handle to a running work-stealing pool. Created by [`Pool::start`] and
/// consumed by [`Pool::join_all`].
pub struct Pool {
    handles: Vec<thread::JoinHandle<Result<()>>>,
}

impl Pool {
    /// Spawn `worker_threads` pool workers.
    ///
    /// `queues.len()` must equal `invokers.len() + 1`: one queue feeds the
    /// first stage, each subsequent queue carries a stage's output into the
    /// next stage (or the sink).
    ///
    /// # Errors
    /// Returns an error if a thread fails to spawn.
    ///
    /// # Panics
    /// Panics if `queues.len() != invokers.len() + 1`.
    #[tracing::instrument(
        name = "pool_spawn_workers",
        skip_all,
        fields(num_stages = invokers.len(), worker_threads = worker_threads),
    )]
    pub fn start(
        invokers: Vec<StageInvoker>,
        queues: Vec<Arc<ErasedQueue>>,
        cancel: &CancelToken,
        worker_threads: usize,
        _global_memory: Arc<MemoryTracker>,
        pool_stats: &Arc<PipelineStats>,
    ) -> Result<Self> {
        assert_eq!(
            queues.len(),
            invokers.len() + 1,
            "queues.len() must be invokers.len() + 1 (one queue per stage pair)",
        );
        assert!(worker_threads > 0, "worker_threads must be >= 1");

        let num_stages = invokers.len();
        let invokers = Arc::new(invokers);
        let queues = Arc::new(queues);

        // One "producers remaining" counter per stage, initialized to
        // `worker_threads`. When a worker observes stage N's input drained
        // and has no held item for stage N, it decrements this counter
        // exactly once. When the counter reaches 0, the worker closes the
        // corresponding output queue (queue[N+1]) so downstream consumers
        // (the next stage or the sink) can observe drain.
        let producers_remaining: Arc<Vec<AtomicUsize>> =
            Arc::new((0..num_stages).map(|_| AtomicUsize::new(worker_threads)).collect());

        let kinds = parallelism_snapshot(&invokers);
        let owned_per_worker = assign_sequential_ownership(worker_threads, &kinds);
        let kinds = Arc::new(kinds);

        let mut handles = Vec::with_capacity(worker_threads);

        for (id, owned) in owned_per_worker.into_iter().enumerate() {
            let invokers_for_worker = invokers.clone();
            let queues_for_worker = queues.clone();
            let cancel_for_worker = cancel.clone();
            let producers_for_worker = producers_remaining.clone();
            let role = if id == 0 && worker_threads > 1 {
                WorkerRole::HeadAffinity
            } else if id + 1 == worker_threads && worker_threads > 1 {
                WorkerRole::TailAffinity
            } else {
                WorkerRole::Generic
            };
            let kinds_for_worker = kinds.clone();
            let stats_for_worker = pool_stats.clone();

            let thread_name = format!("pool_worker_{id}");
            let num_stages = invokers.len();
            let handle =
                thread::Builder::new().name(thread_name.clone()).spawn(move || -> Result<()> {
                    test_support::register_thread(&thread_name);

                    // `stage_exited` lives outside the worker state so that a
                    // panic during `run_worker_loop` does not drop it before
                    // the cleanup guard fires. The guard's `Drop` ensures
                    // every not-yet-exited stage decrements its
                    // `producers_remaining` counter and cascades `close()`
                    // to the output queue if this worker was the last
                    // producer. This fires on every exit path, including
                    // panic unwind.
                    let stage_exited: Arc<Mutex<Vec<bool>>> =
                        Arc::new(Mutex::new(vec![false; num_stages]));

                    let _guard = ExitGuard {
                        stage_exited: stage_exited.clone(),
                        producers_remaining: producers_for_worker.clone(),
                        queues: queues_for_worker.clone(),
                    };

                    let cancel_for_catch = cancel_for_worker.clone();
                    let stage_exited_for_worker = stage_exited.clone();
                    run_and_catch_panic(&cancel_for_catch, move || {
                        let scheduler = Scheduler::new(num_stages, role, owned);
                        let mut state = WorkerState::new(scheduler, &invokers_for_worker);
                        run_worker_loop(
                            &mut state,
                            &invokers_for_worker,
                            &queues_for_worker,
                            &producers_for_worker,
                            &stage_exited_for_worker,
                            &cancel_for_worker,
                            &kinds_for_worker,
                            &stats_for_worker,
                            id,
                        )
                    })
                })?;
            handles.push(handle);
        }

        Ok(Pool { handles })
    }

    /// Join every worker, returning the first error observed. On any error or
    /// panic, cancels the pipeline so sibling workers exit promptly.
    ///
    /// # Errors
    /// Returns the first worker error, or a synthesized error if any worker
    /// panicked.
    pub fn join_all(self, cancel: &CancelToken) -> Result<()> {
        let mut first_err: Option<anyhow::Error> = None;
        for h in self.handles {
            match h.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => {
                    cancel.cancel();
                    if first_err.is_none() {
                        first_err = Some(e);
                    }
                }
                Err(payload) => {
                    cancel.cancel();
                    if first_err.is_none() {
                        // This path only fires if the `run_and_catch_panic`
                        // wrapper itself panicked (which it shouldn't) —
                        // stage-internal panics are already caught by the
                        // per-stage `catch_unwind` in `try_run_stage`.
                        // Extract a readable message instead of `{:?}`-
                        // printing the opaque `Box<dyn Any>`.
                        let msg = extract_panic_message(payload.as_ref());
                        first_err = Some(anyhow::anyhow!("pool worker panicked: {msg}"));
                    }
                }
            }
        }
        if let Some(e) = first_err { Err(e) } else { Ok(()) }
    }
}

/// Try to mark stage `stage_idx` as exited for this worker, decrementing
/// `producers_remaining[stage_idx]` and closing the output queue if this
/// worker was the last producer. Only acts when the stage's input is drained
/// and this worker holds no item for that stage.
fn try_exit_stage(
    stage_idx: usize,
    state: &WorkerState,
    queues: &[Arc<ErasedQueue>],
    producers_remaining: &[AtomicUsize],
    stage_exited: &Mutex<Vec<bool>>,
    exited_snapshot: &[bool],
) {
    if exited_snapshot[stage_idx] {
        return;
    }
    if queues[stage_idx].is_drained() && !state.holds_for_stage(stage_idx) {
        let mut guard = stage_exited.lock().unwrap();
        if !guard[stage_idx] {
            guard[stage_idx] = true;
            let prev = producers_remaining[stage_idx].fetch_sub(1, Ordering::SeqCst);
            if prev == 1 {
                queues[stage_idx + 1].close();
            }
        }
    }
}

/// Main per-worker loop. Runs until the worker has exited every stage it can
/// contribute to (input drained, no held items) or until cancelled.
///
/// Exit-cascade: when a worker observes stage N drained with no held item
/// for stage N, it decrements `producers_remaining[N]`. The last worker to
/// observe drain closes queue[N+1] so the next stage's workers (and the
/// sink) can also observe drain.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn run_worker_loop(
    state: &mut WorkerState,
    invokers: &[StageInvoker],
    queues: &[Arc<ErasedQueue>],
    producers_remaining: &[AtomicUsize],
    stage_exited: &Mutex<Vec<bool>>,
    cancel: &CancelToken,
    parallelism: &[Parallelism],
    pool_stats: &PipelineStats,
    worker_id: usize,
) -> Result<()> {
    // BackpressureState needs a MemoryTracker for its `memory_high` /
    // `memory_drained` flags. Those flags only control drain-mode priority;
    // for now we pass a never-full tracker and let per-queue backpressure
    // drive the scheduler. (The real global tracker lives one level up in
    // the driver — plumbing it through here is a follow-up.)
    let dummy_tracker = Arc::new(MemoryTracker::new(usize::MAX));

    let mut iter_count: u64 = 0;
    let mut cached_bp: Option<BackpressureState> = None;
    let mut last_had_progress = true;

    loop {
        if cancel.is_cancelled() {
            return Ok(());
        }

        // Refresh backpressure sample every 4 iterations, on the first
        // iteration, or whenever we transitioned from making progress to
        // not making progress (so the next snooze is informed by current
        // state rather than stale data).
        let needs_refresh =
            cached_bp.is_none() || iter_count.is_multiple_of(4) || !last_had_progress;
        if needs_refresh {
            let summaries: Vec<QueueSummary> =
                queues.iter().enumerate().map(|(i, q)| q.summarize(i)).collect();
            cached_bp = Some(BackpressureState::new(summaries, &dummy_tracker));
        }
        let bp = cached_bp.as_ref().expect("cached_bp refreshed above");

        // Eager exclusive-step execution (mirrors unified_pipeline's
        // base.rs:4509-4540). Owner workers run their Sequential stages
        // FIRST, before consulting the scheduler's priority list. Without
        // this, the Sequential stage may sit late in the priority list
        // while the owner spends iterations on Parallel stages, starving
        // the Sequential queue and deadlocking the pipeline once the
        // source closes its output.
        let mut any_progress = false;
        let exited_snapshot_eager: Vec<bool> = stage_exited.lock().unwrap().clone();
        let owned: Vec<usize> = state.scheduler.owned_stages().to_vec();
        for stage_idx in owned {
            if exited_snapshot_eager[stage_idx] {
                continue;
            }
            let input = &queues[stage_idx];
            let output = &queues[stage_idx + 1];
            let t0 = Instant::now();
            let result = state.try_run_stage(
                stage_idx,
                &invokers[stage_idx],
                input,
                output,
                cancel,
                pool_stats,
            );
            let elapsed = t0.elapsed();
            match &result {
                StageAttemptResult::Progress => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(stage_idx, StageAttemptOutcome::Success, bp);
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_step_progress(worker_id, stage_idx, elapsed);
                    pool_stats.record_eager_progress(worker_id);
                    state.backoff.reset();
                    any_progress = true;
                }
                StageAttemptResult::BackpressureInput => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(
                        stage_idx,
                        StageAttemptOutcome::BackpressureInput,
                        bp,
                    );
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_bp_input(stage_idx);
                    pool_stats.record_queue_empty_sample(stage_idx);
                }
                StageAttemptResult::BackpressureOutput => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(
                        stage_idx,
                        StageAttemptOutcome::BackpressureOutput,
                        bp,
                    );
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_bp_output(stage_idx);
                }
                StageAttemptResult::SequentialBusy => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(
                        stage_idx,
                        StageAttemptOutcome::SequentialBusy,
                        bp,
                    );
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_seq_busy(stage_idx);
                }
                StageAttemptResult::Drained => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(stage_idx, StageAttemptOutcome::Drained, bp);
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_drained(stage_idx);
                }
                StageAttemptResult::Cancelled | StageAttemptResult::Error(_) => {}
            }
            match result {
                StageAttemptResult::Drained => {
                    try_exit_stage(
                        stage_idx,
                        state,
                        queues,
                        producers_remaining,
                        stage_exited,
                        &exited_snapshot_eager,
                    );
                }
                StageAttemptResult::Cancelled => return Ok(()),
                StageAttemptResult::Error(e) => {
                    cancel.cancel();
                    return Err(e);
                }
                _ => {}
            }
        }

        // Ask the scheduler for stage priorities (ownership-filtered).
        let priorities = state.scheduler.priorities_filtered(bp, parallelism);

        // Snapshot stage_exited for this iteration under the lock; we
        // only update it on observed drain events below.
        let exited_snapshot: Vec<bool> = stage_exited.lock().unwrap().clone();

        for stage_idx in priorities {
            // Skip stages this worker has already exited; those no longer
            // participate in the scheduler for this worker.
            if exited_snapshot[stage_idx] {
                continue;
            }

            let input = &queues[stage_idx];
            let output = &queues[stage_idx + 1];
            let t0 = Instant::now();
            let result = state.try_run_stage(
                stage_idx,
                &invokers[stage_idx],
                input,
                output,
                cancel,
                pool_stats,
            );
            let elapsed = t0.elapsed();
            match &result {
                StageAttemptResult::Progress => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(stage_idx, StageAttemptOutcome::Success, bp);
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_step_progress(worker_id, stage_idx, elapsed);
                    state.backoff.reset();
                    any_progress = true;
                }
                StageAttemptResult::BackpressureInput => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(
                        stage_idx,
                        StageAttemptOutcome::BackpressureInput,
                        bp,
                    );
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_bp_input(stage_idx);
                    pool_stats.record_queue_empty_sample(stage_idx);
                }
                StageAttemptResult::BackpressureOutput => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(
                        stage_idx,
                        StageAttemptOutcome::BackpressureOutput,
                        bp,
                    );
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_bp_output(stage_idx);
                }
                StageAttemptResult::SequentialBusy => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(
                        stage_idx,
                        StageAttemptOutcome::SequentialBusy,
                        bp,
                    );
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_seq_busy(stage_idx);
                }
                StageAttemptResult::Drained => {
                    let was_draining = state.scheduler.is_draining();
                    state.scheduler.record_outcome(stage_idx, StageAttemptOutcome::Drained, bp);
                    if !was_draining && state.scheduler.is_draining() {
                        pool_stats.record_drain_activation();
                    }
                    pool_stats.record_drained(stage_idx);
                }
                StageAttemptResult::Cancelled | StageAttemptResult::Error(_) => {}
            }
            match result {
                StageAttemptResult::Progress => break,
                StageAttemptResult::BackpressureInput
                | StageAttemptResult::BackpressureOutput
                | StageAttemptResult::SequentialBusy => {}
                StageAttemptResult::Drained => {
                    // Input is drained. If we also hold no item for this
                    // stage, this worker has exited the stage. Update the
                    // shared exit flag under the lock so the pool's cleanup
                    // guard sees it, and decrement the atomic producers
                    // counter exactly once. When the counter reaches 0,
                    // close the output queue.
                    try_exit_stage(
                        stage_idx,
                        state,
                        queues,
                        producers_remaining,
                        stage_exited,
                        &exited_snapshot,
                    );
                }
                StageAttemptResult::Cancelled => return Ok(()),
                StageAttemptResult::Error(e) => {
                    cancel.cancel();
                    return Err(e);
                }
            }
        }

        // Drain-check for Sequential stages this worker does not own.
        // Non-owner workers are still registered as producers in
        // `producers_remaining` (initialized to `worker_threads`). When the
        // owner finishes and the input queue drains, non-owners must also
        // decrement their counter so the output queue eventually closes.
        // They never attempt the stage, so this sweep is their only exit path.
        for (stage_idx, kind) in parallelism.iter().enumerate() {
            if *kind == Parallelism::Sequential && !state.scheduler.owns(stage_idx) {
                try_exit_stage(
                    stage_idx,
                    state,
                    queues,
                    producers_remaining,
                    stage_exited,
                    &exited_snapshot,
                );
            }
        }

        last_had_progress = any_progress;
        iter_count = iter_count.wrapping_add(1);

        if any_progress {
            continue;
        }

        // Exit condition: this worker has exited every stage AND has no
        // held items outstanding.
        let all_exited = stage_exited.lock().unwrap().iter().all(|d| *d);
        if all_exited && !state.has_held_items() {
            return Ok(());
        }

        // No progress this iteration — back off to avoid burning a core.
        let snooze_start = Instant::now();
        state.backoff.snooze();
        pool_stats.record_worker_idle(worker_id, snooze_start.elapsed());
    }
}

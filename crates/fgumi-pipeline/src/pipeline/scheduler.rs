//! Sticky-consensus scheduler for Zone 3 pipeline stages.
//!
//! Most workers stay on Consensus (the CPU bottleneck). One dedicated writer prevents Write
//! contention. When a worker's current step has no work, it chases downstream to drain the
//! pipeline.
//!
//! Each worker owns a [`WorkerState`] with held-item slots, enabling non-blocking push:
//! when an output queue is full, the item is stored in the worker's state and retried on the
//! next iteration. This eliminates the per-item spin-wait of the old blocking push.

use std::io::Write;
use std::panic::AssertUnwindSafe;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use anyhow::{Result, anyhow};
use crossbeam_channel::{Receiver, Sender};

use crate::pipeline::WorkQueue;
use crate::pipeline::reorder::ReorderBuffer;
use crate::pipeline::stats::StageStats;
use crate::pipeline::worker_util::{WorkerBackoff, panic_message};
use crate::stage::{PipelineStage, SequencedItem};
use crate::stages::compress::CompressedBatch;
use crate::stages::group_assign::{MiGroup, PositionGroupBatch};

// ============================================================================
// Constants
// ============================================================================

/// Number of pipeline stages.
const NUM_STAGES: usize = 5;
/// Number of stages that can hold items (all except Write).
const NUM_HOLDABLE: usize = 4;

/// Stage indices.
const GROUP_ASSIGN: usize = 0;
const CONSENSUS: usize = 1;
const FILTER: usize = 2;
const COMPRESS: usize = 3;
const WRITE: usize = 4;

// ============================================================================
// StepResult
// ============================================================================

/// Result of attempting a pipeline step.
#[derive(Debug)]
pub enum StepResult {
    /// Successfully processed and pushed an item.
    Success,
    /// Output queue is full. Item is held in worker state.
    OutputFull,
    /// Input queue is empty or closed.
    InputEmpty,
    /// Processing error. The caller should propagate this and cancel.
    Error(anyhow::Error),
}

// ============================================================================
// PipelineStep / BackpressureState
// ============================================================================

/// The five pipeline stages in the runall pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipelineStep {
    GroupAssign = 0,
    Consensus = 1,
    Filter = 2,
    Compress = 3,
    Write = 4,
}

impl PipelineStep {
    /// Convert an index (0..5) to a `PipelineStep`.
    #[cfg(test)]
    fn from_index(idx: usize) -> Self {
        match idx {
            0 => Self::GroupAssign,
            1 => Self::Consensus,
            2 => Self::Filter,
            3 => Self::Compress,
            4 => Self::Write,
            _ => panic!("invalid PipelineStep index: {idx}"),
        }
    }
}

/// Backpressure signals sampled from queue state each iteration.
pub struct BackpressureState {
    /// True when `q_compressed` (Compress output) exceeds high-water mark.
    pub output_high: bool,
    /// True when `q_mi_batch` memory tracker is at or above limit.
    pub memory_high: bool,
}

impl BackpressureState {
    /// Sample current backpressure state from the `StageGraph`'s queues.
    #[must_use]
    pub fn sample(graph: &StageGraph) -> Self {
        Self {
            output_high: graph.q_compressed.fill_ratio() > super::BACKPRESSURE_HIGH_WATERMARK,
            memory_high: graph.q_mi_batch.memory_ratio() >= 1.0,
        }
    }
}

// ============================================================================
// PipelineScheduler — sticky-consensus strategy
// ============================================================================

/// Thread role within the runall pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkerRole {
    /// Worker N-1: tries Write first, then Consensus.
    WritePreferred,
    /// Worker 1: tries `GroupAssign` first, then Consensus (only when >= 4 workers).
    GroupAssignPreferred,
    /// Workers 2..N-2: sticky on Consensus, chase downstream on idle.
    ConsensusWorker,
}

/// Per-worker scheduler for the runall pipeline.
///
/// Uses a sticky-consensus strategy: most workers stay on Consensus (the CPU bottleneck,
/// 49% of time). One dedicated writer prevents Write contention. When a worker's current
/// step has no work, it chases downstream to drain the pipeline.
pub struct PipelineScheduler {
    /// Current preferred step — stays on last successful step.
    current_step: PipelineStep,
    /// Priority output buffer.
    priority_buffer: [PipelineStep; NUM_STAGES],
    /// Thread role.
    pub role: WorkerRole,
}

impl PipelineScheduler {
    /// Create a new scheduler for the given worker.
    #[must_use]
    pub fn new(worker_id: usize, num_workers: usize) -> Self {
        let role = if worker_id == num_workers.saturating_sub(1) {
            WorkerRole::WritePreferred
        } else if worker_id == 1 && num_workers >= 4 {
            WorkerRole::GroupAssignPreferred
        } else {
            WorkerRole::ConsensusWorker
        };

        let current_step = match role {
            WorkerRole::WritePreferred => PipelineStep::Compress,
            WorkerRole::GroupAssignPreferred => PipelineStep::GroupAssign,
            WorkerRole::ConsensusWorker => PipelineStep::Consensus,
        };

        Self { current_step, priority_buffer: [PipelineStep::GroupAssign; NUM_STAGES], role }
    }

    /// Build priority list based on current state and backpressure signals.
    ///
    /// Priority logic:
    /// 1. Exclusive role step first (if any)
    /// 2. Current step (sticky on last success)
    /// 3. Chase downstream: Filter -> Compress -> Write
    /// 4. Chase upstream: `GroupAssign`
    /// 5. Remaining steps
    ///
    /// Backpressure override: if `output_high`, move Compress to front (after exclusive role).
    pub fn get_priorities(&mut self, bp: &BackpressureState) -> &[PipelineStep; NUM_STAGES] {
        let mut buf: [PipelineStep; NUM_STAGES] = [PipelineStep::GroupAssign; NUM_STAGES];
        let mut len = 0;

        let mut push = |step: PipelineStep| {
            if !buf[..len].contains(&step) {
                buf[len] = step;
                len += 1;
            }
        };

        // 1. Exclusive role step first.
        match self.role {
            WorkerRole::WritePreferred => push(PipelineStep::Write),
            WorkerRole::GroupAssignPreferred => push(PipelineStep::GroupAssign),
            WorkerRole::ConsensusWorker => {}
        }

        // Backpressure override: if output_high, Compress goes right after exclusive role.
        if bp.output_high {
            push(PipelineStep::Compress);
        }

        // 2. Current step (sticky).
        push(self.current_step);

        // 3. Chase downstream: Filter -> Compress -> Write.
        push(PipelineStep::Filter);
        push(PipelineStep::Compress);
        push(PipelineStep::Write);

        // 4. Chase upstream: Consensus -> GroupAssign.
        push(PipelineStep::Consensus);
        push(PipelineStep::GroupAssign);

        self.priority_buffer = buf;
        &self.priority_buffer
    }

    /// Record outcome of attempting a step.
    ///
    /// On success, the current step becomes sticky on the successful step.
    /// Special pivots: `WritePreferred` pivots to Compress after Write,
    /// `GroupAssignPreferred` pivots to Consensus after `GroupAssign`.
    pub fn record_outcome(&mut self, step: PipelineStep, success: bool) {
        if success {
            self.current_step = match self.role {
                WorkerRole::WritePreferred if step == PipelineStep::Write => PipelineStep::Compress,
                WorkerRole::GroupAssignPreferred if step == PipelineStep::GroupAssign => {
                    PipelineStep::Consensus
                }
                _ => step,
            };
        }
        // On failure: leave current_step unchanged; the priority list naturally tries
        // the next step.
    }
}

// ============================================================================
// WorkerStats — per-worker comprehensive statistics
// ============================================================================

/// Per-worker statistics for pipeline diagnostics.
///
/// Each worker owns its own `WorkerStats` instance — no cross-thread
/// contention on stat counters.
pub struct WorkerStats {
    /// Per-stage: number of times `try_step_*` was called.
    pub stage_attempts: [AtomicU64; NUM_STAGES],
    /// Per-stage: number of successful `try_step_*` calls.
    pub stage_successes: [AtomicU64; NUM_STAGES],
    /// Per-stage: cumulative nanoseconds spent in `process()` calls.
    pub stage_nanos: [AtomicU64; NUM_STAGES],
    /// Per-stage: cumulative items processed.
    pub stage_items: [AtomicU64; NUM_STAGES],

    /// Number of failed `try_lock` attempts on the write mutex.
    pub write_lock_contention: AtomicU64,
    /// Per-stage: number of times `push()` returned full.
    pub push_full_count: [AtomicU64; NUM_STAGES],
    /// Per-stage: number of times `try_pop()` returned None.
    pub pop_empty_count: [AtomicU64; NUM_STAGES],

    /// Per-holdable-stage: number of times an item was held.
    pub held_count: [AtomicU64; NUM_HOLDABLE],

    /// Number of full iterations where no work was found.
    pub idle_cycles: AtomicU64,
    /// Cumulative nanoseconds spent in backoff sleep.
    pub idle_nanos: AtomicU64,
}

impl WorkerStats {
    /// Create a new zeroed stats instance.
    fn new() -> Self {
        Self {
            stage_attempts: std::array::from_fn(|_| AtomicU64::new(0)),
            stage_successes: std::array::from_fn(|_| AtomicU64::new(0)),
            stage_nanos: std::array::from_fn(|_| AtomicU64::new(0)),
            stage_items: std::array::from_fn(|_| AtomicU64::new(0)),
            write_lock_contention: AtomicU64::new(0),
            push_full_count: std::array::from_fn(|_| AtomicU64::new(0)),
            pop_empty_count: std::array::from_fn(|_| AtomicU64::new(0)),
            held_count: std::array::from_fn(|_| AtomicU64::new(0)),
            idle_cycles: AtomicU64::new(0),
            idle_nanos: AtomicU64::new(0),
        }
    }
}

// ============================================================================
// WorkerState — per-worker mutable state
// ============================================================================

/// Per-worker mutable state for all stages.
///
/// Each worker thread owns one `WorkerState`. Held items are stored here
/// when an output queue is full, enabling non-blocking push semantics.
pub struct WorkerState {
    /// Worker ID (0-indexed).
    pub worker_id: usize,

    // Held output items for each inter-stage queue. Each held item is paired with
    // the input memory estimate needed to call `release_memory()` on the input queue
    // when the held item is finally pushed.
    held_group_assign: Option<(SequencedItem<Vec<MiGroup>>, usize)>,
    held_consensus: Option<(SequencedItem<fgumi_consensus::caller::ConsensusOutput>, usize)>,
    held_filter: Option<(SequencedItem<Vec<u8>>, usize)>,
    held_compress: Option<(SequencedItem<CompressedBatch>, usize)>,

    /// Adaptive backoff with jitter.
    backoff: WorkerBackoff,

    /// Sticky-consensus scheduler — determines priority order each iteration.
    pub scheduler: PipelineScheduler,
    /// Per-worker stats.
    pub stats: WorkerStats,
}

impl WorkerState {
    /// Create a new `WorkerState` for the given worker.
    #[must_use]
    pub fn new(worker_id: usize, num_workers: usize) -> Self {
        Self {
            worker_id,
            held_group_assign: None,
            held_consensus: None,
            held_filter: None,
            held_compress: None,
            backoff: WorkerBackoff::new(worker_id),
            scheduler: PipelineScheduler::new(worker_id, num_workers),
            stats: WorkerStats::new(),
        }
    }

    /// Returns `true` if this worker is holding any items.
    pub fn has_any_held_items(&self) -> bool {
        self.held_group_assign.is_some()
            || self.held_consensus.is_some()
            || self.held_filter.is_some()
            || self.held_compress.is_some()
    }

    /// Reset backoff to zero on successful work.
    pub fn reset_backoff(&mut self) {
        self.backoff.reset();
    }

    /// Sleep for the current backoff duration with per-worker jitter.
    pub fn sleep_backoff(&mut self) {
        self.backoff.sleep();
    }

    /// Increase backoff exponentially.
    pub fn increase_backoff(&mut self) {
        self.backoff.increase();
    }

    /// Log comprehensive stats for this worker.
    #[expect(clippy::cast_precision_loss, reason = "nanos fits in f64 for display")]
    pub fn log_stats(&self) {
        let role = match self.scheduler.role {
            WorkerRole::WritePreferred => " (write-preferred)",
            WorkerRole::GroupAssignPreferred => " (group-assign-preferred)",
            WorkerRole::ConsensusWorker => " (consensus-worker)",
        };
        log::info!("  Worker {}{role}:", self.worker_id);

        let stage_names = ["group_assign", "consensus", "filter", "compress", "write"];
        for (idx, name) in stage_names.iter().enumerate() {
            let attempts = self.stats.stage_attempts[idx].load(Ordering::Relaxed);
            let successes = self.stats.stage_successes[idx].load(Ordering::Relaxed);
            let nanos = self.stats.stage_nanos[idx].load(Ordering::Relaxed);
            let items = self.stats.stage_items[idx].load(Ordering::Relaxed);
            let empty = self.stats.pop_empty_count[idx].load(Ordering::Relaxed);
            let full = self.stats.push_full_count[idx].load(Ordering::Relaxed);
            let secs = nanos as f64 / 1_000_000_000.0;
            log::info!(
                "    {name:<14} {attempts:>8} attempts, {successes:>8} success, \
                 {items:>8} items, {empty:>6} empty, {full:>6} full, {secs:.2}s"
            );
        }

        let contention = self.stats.write_lock_contention.load(Ordering::Relaxed);
        if contention > 0 {
            log::info!("    write_contention: {contention}");
        }

        let held_names = ["group_assign", "consensus", "filter", "compress"];
        for (idx, name) in held_names.iter().enumerate() {
            let held = self.stats.held_count[idx].load(Ordering::Relaxed);
            if held > 0 {
                log::info!("    {name:<14} {held:>8} held");
            }
        }

        let idle_cycles = self.stats.idle_cycles.load(Ordering::Relaxed);
        let idle_nanos = self.stats.idle_nanos.load(Ordering::Relaxed);
        let idle_secs = idle_nanos as f64 / 1_000_000_000.0;
        log::info!("    idle: {idle_cycles} cycles, {idle_secs:.3}s");
    }
}

// ============================================================================
// WriteState — shared mutable state for the write stage
// ============================================================================

/// Shared mutable state for the write stage, protected by a [`Mutex`].
///
/// Any worker thread (or the calling thread) can acquire the lock and perform writes.
pub struct WriteState {
    /// Output file handle.
    pub(crate) file: std::fs::File,
    /// Reorder buffer ensuring compressed batches are written in sequence order.
    pub(crate) reorder_buf: ReorderBuffer<CompressedBatch>,
    /// Counter tracking how many compressed batches have been written.
    pub(crate) batches_written: u64,
}

// ============================================================================
// StageGraph — typed stage graph replacing opaque closures
// ============================================================================

/// The typed stage graph for the runall pipeline.
///
/// Holds `Arc` references to all stages and queues. Shared (immutable) across
/// all worker threads. Workers iterate stages via the `try_step_*` functions,
/// using their own [`WorkerState`] for held items.
pub struct StageGraph {
    // Stages.
    pub(crate) group_assign:
        Arc<dyn PipelineStage<Input = PositionGroupBatch, Output = Vec<MiGroup>>>,
    pub(crate) consensus: Arc<
        dyn PipelineStage<Input = Vec<MiGroup>, Output = fgumi_consensus::caller::ConsensusOutput>,
    >,
    pub(crate) filter:
        Arc<dyn PipelineStage<Input = fgumi_consensus::caller::ConsensusOutput, Output = Vec<u8>>>,
    pub(crate) compress: Arc<dyn PipelineStage<Input = Vec<u8>, Output = CompressedBatch>>,

    // Queues (ordered upstream to downstream).
    pub(crate) q_position_batch: WorkQueue<SequencedItem<PositionGroupBatch>>,
    pub(crate) q_mi_batch: WorkQueue<SequencedItem<Vec<MiGroup>>>,
    pub(crate) q_consensus: WorkQueue<SequencedItem<fgumi_consensus::caller::ConsensusOutput>>,
    pub(crate) q_filtered: WorkQueue<SequencedItem<Vec<u8>>>,
    pub(crate) q_compressed: WorkQueue<SequencedItem<CompressedBatch>>,

    // Write state (mutex-protected, highest priority).
    pub(crate) write_state: Arc<Mutex<WriteState>>,

    // Per-stage stats.
    pub(crate) stats: [Arc<StageStats>; NUM_STAGES],

    // Cancellation.
    pub(crate) cancel: Arc<AtomicBool>,
}

// ============================================================================
// try_push_held_* — advance held items from previous iteration
// ============================================================================

/// Try to push a held `group_assign` output. Returns true if pushed.
pub(crate) fn try_push_held_group_assign(graph: &StageGraph, worker: &mut WorkerState) -> bool {
    let Some((held, input_mem)) = worker.held_group_assign.take() else {
        return false;
    };
    let output_mem = held.memory_estimate;
    if let Some(returned) = graph.q_mi_batch.push(held, output_mem) {
        // Still full — put it back.
        worker.held_group_assign = Some((returned, input_mem));
        false
    } else {
        // Successfully pushed.
        graph.q_position_batch.complete_in_flight();
        graph.q_position_batch.release_memory(input_mem);
        true
    }
}

/// Try to push a held consensus output. Returns true if pushed.
pub(crate) fn try_push_held_consensus(graph: &StageGraph, worker: &mut WorkerState) -> bool {
    let Some((held, input_mem)) = worker.held_consensus.take() else {
        return false;
    };
    let output_mem = held.memory_estimate;
    if let Some(returned) = graph.q_consensus.push(held, output_mem) {
        worker.held_consensus = Some((returned, input_mem));
        false
    } else {
        graph.q_mi_batch.complete_in_flight();
        graph.q_mi_batch.release_memory(input_mem);
        true
    }
}

/// Try to push a held filter output. Returns true if pushed.
pub(crate) fn try_push_held_filter(graph: &StageGraph, worker: &mut WorkerState) -> bool {
    let Some((held, input_mem)) = worker.held_filter.take() else {
        return false;
    };
    let output_mem = held.memory_estimate;
    if let Some(returned) = graph.q_filtered.push(held, output_mem) {
        worker.held_filter = Some((returned, input_mem));
        false
    } else {
        graph.q_consensus.complete_in_flight();
        graph.q_consensus.release_memory(input_mem);
        true
    }
}

/// Try to push a held compress output. Returns true if pushed.
pub(crate) fn try_push_held_compress(graph: &StageGraph, worker: &mut WorkerState) -> bool {
    let Some((held, input_mem)) = worker.held_compress.take() else {
        return false;
    };
    let output_mem = held.memory_estimate;
    if let Some(returned) = graph.q_compressed.push(held, output_mem) {
        worker.held_compress = Some((returned, input_mem));
        false
    } else {
        graph.q_filtered.complete_in_flight();
        graph.q_filtered.release_memory(input_mem);
        true
    }
}

// ============================================================================
// try_step_* — one per pipeline stage
// ============================================================================

/// Attempt one unit of group-assign work.
pub(crate) fn try_step_group_assign(graph: &StageGraph, worker: &mut WorkerState) -> StepResult {
    worker.stats.stage_attempts[GROUP_ASSIGN].fetch_add(1, Ordering::Relaxed);

    // 1. If holding an item from a previous iteration, try to push it first.
    if let Some((held, input_mem)) = worker.held_group_assign.take() {
        let output_mem = held.memory_estimate;
        if let Some(returned) = graph.q_mi_batch.push(held, output_mem) {
            worker.held_group_assign = Some((returned, input_mem));
            worker.stats.push_full_count[GROUP_ASSIGN].fetch_add(1, Ordering::Relaxed);
            return StepResult::OutputFull;
        }
        graph.q_position_batch.complete_in_flight();
        graph.q_position_batch.release_memory(input_mem);
        // Held item pushed — fall through to try new input.
    }

    // 2. Pop input.
    let Some(input) = graph.q_position_batch.try_pop() else {
        worker.stats.pop_empty_count[GROUP_ASSIGN].fetch_add(1, Ordering::Relaxed);
        if graph.q_position_batch.is_closed_and_empty() {
            graph.q_mi_batch.close();
        }
        return StepResult::InputEmpty;
    };
    let input_mem = input.memory_estimate;

    // 3. Process.
    let start = Instant::now();
    let result = match graph.group_assign.process(input.item) {
        Ok(r) => r,
        Err(e) => {
            graph.q_position_batch.complete_in_flight();
            graph.q_position_batch.release_memory(input_mem);
            return StepResult::Error(e);
        }
    };
    #[expect(clippy::cast_possible_truncation, reason = "elapsed nanos won't exceed u64")]
    let elapsed_nanos = start.elapsed().as_nanos() as u64;
    worker.stats.stage_nanos[GROUP_ASSIGN].fetch_add(elapsed_nanos, Ordering::Relaxed);
    graph.stats[GROUP_ASSIGN].record(elapsed_nanos, 1);

    let output_mem = graph.group_assign.output_memory_estimate(&result);
    let output = SequencedItem::new(input.seq, result, output_mem);

    // 4. Push or hold.
    if let Some(returned) = graph.q_mi_batch.push(output, output_mem) {
        worker.held_group_assign = Some((returned, input_mem));
        worker.stats.push_full_count[GROUP_ASSIGN].fetch_add(1, Ordering::Relaxed);
        worker.stats.held_count[GROUP_ASSIGN].fetch_add(1, Ordering::Relaxed);
        StepResult::OutputFull
    } else {
        graph.q_position_batch.complete_in_flight();
        graph.q_position_batch.release_memory(input_mem);
        worker.stats.stage_successes[GROUP_ASSIGN].fetch_add(1, Ordering::Relaxed);
        worker.stats.stage_items[GROUP_ASSIGN].fetch_add(1, Ordering::Relaxed);
        StepResult::Success
    }
}

/// Attempt one unit of consensus work.
pub(crate) fn try_step_consensus(graph: &StageGraph, worker: &mut WorkerState) -> StepResult {
    worker.stats.stage_attempts[CONSENSUS].fetch_add(1, Ordering::Relaxed);

    if let Some((held, input_mem)) = worker.held_consensus.take() {
        let output_mem = held.memory_estimate;
        if let Some(returned) = graph.q_consensus.push(held, output_mem) {
            worker.held_consensus = Some((returned, input_mem));
            worker.stats.push_full_count[CONSENSUS].fetch_add(1, Ordering::Relaxed);
            return StepResult::OutputFull;
        }
        graph.q_mi_batch.complete_in_flight();
        graph.q_mi_batch.release_memory(input_mem);
    }

    let Some(input) = graph.q_mi_batch.try_pop() else {
        worker.stats.pop_empty_count[CONSENSUS].fetch_add(1, Ordering::Relaxed);
        if graph.q_mi_batch.is_closed_and_empty() {
            graph.q_consensus.close();
        }
        return StepResult::InputEmpty;
    };
    let input_mem = input.memory_estimate;

    let start = Instant::now();
    let result = match graph.consensus.process(input.item) {
        Ok(r) => r,
        Err(e) => {
            graph.q_mi_batch.complete_in_flight();
            graph.q_mi_batch.release_memory(input_mem);
            return StepResult::Error(e);
        }
    };
    #[expect(clippy::cast_possible_truncation, reason = "elapsed nanos won't exceed u64")]
    let elapsed_nanos = start.elapsed().as_nanos() as u64;
    worker.stats.stage_nanos[CONSENSUS].fetch_add(elapsed_nanos, Ordering::Relaxed);
    graph.stats[CONSENSUS].record(elapsed_nanos, 1);

    let output_mem = graph.consensus.output_memory_estimate(&result);
    let output = SequencedItem::new(input.seq, result, output_mem);

    if let Some(returned) = graph.q_consensus.push(output, output_mem) {
        worker.held_consensus = Some((returned, input_mem));
        worker.stats.push_full_count[CONSENSUS].fetch_add(1, Ordering::Relaxed);
        worker.stats.held_count[CONSENSUS].fetch_add(1, Ordering::Relaxed);
        StepResult::OutputFull
    } else {
        graph.q_mi_batch.complete_in_flight();
        graph.q_mi_batch.release_memory(input_mem);
        worker.stats.stage_successes[CONSENSUS].fetch_add(1, Ordering::Relaxed);
        worker.stats.stage_items[CONSENSUS].fetch_add(1, Ordering::Relaxed);
        StepResult::Success
    }
}

/// Attempt one unit of filter work.
pub(crate) fn try_step_filter(graph: &StageGraph, worker: &mut WorkerState) -> StepResult {
    worker.stats.stage_attempts[FILTER].fetch_add(1, Ordering::Relaxed);

    if let Some((held, input_mem)) = worker.held_filter.take() {
        let output_mem = held.memory_estimate;
        if let Some(returned) = graph.q_filtered.push(held, output_mem) {
            worker.held_filter = Some((returned, input_mem));
            worker.stats.push_full_count[FILTER].fetch_add(1, Ordering::Relaxed);
            return StepResult::OutputFull;
        }
        graph.q_consensus.complete_in_flight();
        graph.q_consensus.release_memory(input_mem);
    }

    let Some(input) = graph.q_consensus.try_pop() else {
        worker.stats.pop_empty_count[FILTER].fetch_add(1, Ordering::Relaxed);
        if graph.q_consensus.is_closed_and_empty() {
            graph.q_filtered.close();
        }
        return StepResult::InputEmpty;
    };
    let input_mem = input.memory_estimate;

    let start = Instant::now();
    let result = match graph.filter.process(input.item) {
        Ok(r) => r,
        Err(e) => {
            graph.q_consensus.complete_in_flight();
            graph.q_consensus.release_memory(input_mem);
            return StepResult::Error(e);
        }
    };
    #[expect(clippy::cast_possible_truncation, reason = "elapsed nanos won't exceed u64")]
    let elapsed_nanos = start.elapsed().as_nanos() as u64;
    worker.stats.stage_nanos[FILTER].fetch_add(elapsed_nanos, Ordering::Relaxed);
    graph.stats[FILTER].record(elapsed_nanos, 1);

    let output_mem = graph.filter.output_memory_estimate(&result);
    let output = SequencedItem::new(input.seq, result, output_mem);

    if let Some(returned) = graph.q_filtered.push(output, output_mem) {
        worker.held_filter = Some((returned, input_mem));
        worker.stats.push_full_count[FILTER].fetch_add(1, Ordering::Relaxed);
        worker.stats.held_count[FILTER].fetch_add(1, Ordering::Relaxed);
        StepResult::OutputFull
    } else {
        graph.q_consensus.complete_in_flight();
        graph.q_consensus.release_memory(input_mem);
        worker.stats.stage_successes[FILTER].fetch_add(1, Ordering::Relaxed);
        worker.stats.stage_items[FILTER].fetch_add(1, Ordering::Relaxed);
        StepResult::Success
    }
}

/// Attempt one unit of compress work.
pub(crate) fn try_step_compress(graph: &StageGraph, worker: &mut WorkerState) -> StepResult {
    worker.stats.stage_attempts[COMPRESS].fetch_add(1, Ordering::Relaxed);

    if let Some((held, input_mem)) = worker.held_compress.take() {
        let output_mem = held.memory_estimate;
        if let Some(returned) = graph.q_compressed.push(held, output_mem) {
            worker.held_compress = Some((returned, input_mem));
            worker.stats.push_full_count[COMPRESS].fetch_add(1, Ordering::Relaxed);
            return StepResult::OutputFull;
        }
        graph.q_filtered.complete_in_flight();
        graph.q_filtered.release_memory(input_mem);
    }

    let Some(input) = graph.q_filtered.try_pop() else {
        worker.stats.pop_empty_count[COMPRESS].fetch_add(1, Ordering::Relaxed);
        if graph.q_filtered.is_closed_and_empty() {
            graph.q_compressed.close();
        }
        return StepResult::InputEmpty;
    };
    let input_mem = input.memory_estimate;

    let start = Instant::now();
    let result = match graph.compress.process(input.item) {
        Ok(r) => r,
        Err(e) => {
            graph.q_filtered.complete_in_flight();
            graph.q_filtered.release_memory(input_mem);
            return StepResult::Error(e);
        }
    };
    #[expect(clippy::cast_possible_truncation, reason = "elapsed nanos won't exceed u64")]
    let elapsed_nanos = start.elapsed().as_nanos() as u64;
    worker.stats.stage_nanos[COMPRESS].fetch_add(elapsed_nanos, Ordering::Relaxed);
    graph.stats[COMPRESS].record(elapsed_nanos, 1);

    let output_mem = graph.compress.output_memory_estimate(&result);
    let output = SequencedItem::new(input.seq, result, output_mem);

    if let Some(returned) = graph.q_compressed.push(output, output_mem) {
        worker.held_compress = Some((returned, input_mem));
        worker.stats.push_full_count[COMPRESS].fetch_add(1, Ordering::Relaxed);
        worker.stats.held_count[COMPRESS].fetch_add(1, Ordering::Relaxed);
        StepResult::OutputFull
    } else {
        graph.q_filtered.complete_in_flight();
        graph.q_filtered.release_memory(input_mem);
        worker.stats.stage_successes[COMPRESS].fetch_add(1, Ordering::Relaxed);
        worker.stats.stage_items[COMPRESS].fetch_add(1, Ordering::Relaxed);
        StepResult::Success
    }
}

/// Attempt one unit of write work.
///
/// The write stage pops ALL available items from `q_compressed` in a single lock
/// acquisition, inserts them into the reorder buffer, and drains consecutive ready
/// batches to the output file. This batch-drain pattern reduces lock contention.
pub(crate) fn try_step_write(graph: &StageGraph, worker: &mut WorkerState) -> StepResult {
    worker.stats.stage_attempts[WRITE].fetch_add(1, Ordering::Relaxed);

    // Try to acquire write lock non-blocking.
    let mut state = match graph.write_state.try_lock() {
        Ok(guard) => guard,
        Err(std::sync::TryLockError::WouldBlock) => {
            worker.stats.write_lock_contention.fetch_add(1, Ordering::Relaxed);
            return StepResult::InputEmpty;
        }
        Err(std::sync::TryLockError::Poisoned(e)) => {
            return StepResult::Error(anyhow!("write state lock poisoned: {e}"));
        }
    };

    // Drain ALL available items from the compressed queue.
    let mut popped_any = false;
    while let Some(sequenced) = graph.q_compressed.try_pop() {
        if graph.cancel.load(Ordering::Acquire) {
            graph.q_compressed.complete_in_flight();
            return StepResult::InputEmpty;
        }
        let mem: usize = sequenced.item.blocks.iter().map(|b| b.data.len()).sum();
        state.reorder_buf.insert(sequenced.seq, sequenced.item);
        graph.q_compressed.release_memory(mem);
        graph.q_compressed.complete_in_flight();
        popped_any = true;
    }

    if !popped_any && state.reorder_buf.is_empty() {
        worker.stats.pop_empty_count[WRITE].fetch_add(1, Ordering::Relaxed);
        return StepResult::InputEmpty;
    }

    let start = Instant::now();

    // Drain consecutive ready batches and write them.
    let mut batches_written: u64 = 0;
    for batch in state.reorder_buf.drain_ready() {
        for block in &batch.blocks {
            if let Err(e) = state.file.write_all(&block.data) {
                return StepResult::Error(
                    anyhow::Error::new(e).context("writing compressed BGZF block"),
                );
            }
        }
        state.batches_written += 1;
        batches_written += 1;
    }

    #[expect(clippy::cast_possible_truncation, reason = "elapsed nanos won't exceed u64")]
    let elapsed_nanos = start.elapsed().as_nanos() as u64;
    worker.stats.stage_nanos[WRITE].fetch_add(elapsed_nanos, Ordering::Relaxed);
    graph.stats[WRITE].record(elapsed_nanos, batches_written);

    if popped_any || batches_written > 0 {
        worker.stats.stage_successes[WRITE].fetch_add(1, Ordering::Relaxed);
        worker.stats.stage_items[WRITE].fetch_add(batches_written, Ordering::Relaxed);
        StepResult::Success
    } else {
        StepResult::InputEmpty
    }
}

// ============================================================================
// Worker loop
// ============================================================================

/// The main loop executed by each worker thread.
///
/// Each iteration: (1) try to push held items, (2) sample backpressure and get
/// priority order, (3) execute one stage in priority order, (4) backoff if idle.
fn worker_loop(
    graph: &StageGraph,
    worker: &mut WorkerState,
    cancel: &AtomicBool,
    error_tx: &Sender<anyhow::Error>,
) {
    loop {
        // Only exit when cancelled AND no held items remain.
        if cancel.load(Ordering::Acquire) && !worker.has_any_held_items() {
            break;
        }

        let mut did_work = false;

        // 1. Try to advance any held items first (most downstream first).
        did_work |= try_push_held_compress(graph, worker);
        did_work |= try_push_held_filter(graph, worker);
        did_work |= try_push_held_consensus(graph, worker);
        did_work |= try_push_held_group_assign(graph, worker);

        // 2. Sample backpressure and get priority order.
        let bp = BackpressureState::sample(graph);
        let priorities = *worker.scheduler.get_priorities(&bp);

        // 3. Execute one stage in priority order.
        for &stage in &priorities {
            if cancel.load(Ordering::Acquire) {
                break;
            }

            let result = match stage {
                PipelineStep::Write => try_step_write(graph, worker),
                PipelineStep::Compress => try_step_compress(graph, worker),
                PipelineStep::Filter => try_step_filter(graph, worker),
                PipelineStep::Consensus => try_step_consensus(graph, worker),
                PipelineStep::GroupAssign => try_step_group_assign(graph, worker),
            };

            match result {
                StepResult::Success => {
                    worker.scheduler.record_outcome(stage, true);
                    did_work = true;
                    break;
                }
                StepResult::OutputFull => {
                    worker.scheduler.record_outcome(stage, false);
                }
                StepResult::InputEmpty => {}
                StepResult::Error(e) => {
                    let _ = error_tx.send(e);
                    cancel.store(true, Ordering::Release);
                    return;
                }
            }
        }

        // 4. Adaptive backoff if no work done.
        if did_work {
            worker.reset_backoff();
        } else {
            worker.stats.idle_cycles.fetch_add(1, Ordering::Relaxed);
            let start = Instant::now();
            worker.sleep_backoff();
            #[expect(clippy::cast_possible_truncation, reason = "sleep nanos won't exceed u64")]
            let idle_nanos = start.elapsed().as_nanos() as u64;
            worker.stats.idle_nanos.fetch_add(idle_nanos, Ordering::Relaxed);
            worker.increase_backoff();
        }
    }
}

// ============================================================================
// Scheduler
// ============================================================================

/// Shared worker state slot used to return stats from a worker thread after it completes.
pub type WorkerStateSlot = Arc<Mutex<Option<WorkerState>>>;

/// A work-stealing scheduler that runs Zone 3 pipeline stages across a pool of worker
/// threads using a typed [`StageGraph`] and per-worker [`WorkerState`] with held-item slots.
pub struct Scheduler {
    graph: Arc<StageGraph>,
    cancel: Arc<AtomicBool>,
    error_tx: Sender<anyhow::Error>,
    error_rx: Receiver<anyhow::Error>,
}

impl Scheduler {
    /// Create a new scheduler from a [`StageGraph`].
    #[must_use]
    pub fn new(graph: StageGraph) -> Self {
        let cancel = Arc::clone(&graph.cancel);
        let (error_tx, error_rx) = crossbeam_channel::unbounded();
        Self { graph: Arc::new(graph), cancel, error_tx, error_rx }
    }

    /// Start `num_workers` worker threads. Returns their join handles and
    /// the per-worker `WorkerState` instances (for stats collection after join).
    ///
    /// Each worker runs the held-item priority loop until cancelled and all held items
    /// are flushed.
    ///
    /// # Panics
    ///
    /// Panics if the internal worker-state mutex is poisoned.
    #[must_use]
    pub fn start(&self, num_workers: usize) -> (Vec<thread::JoinHandle<()>>, Vec<WorkerStateSlot>) {
        let mut handles = Vec::with_capacity(num_workers);
        let mut worker_states = Vec::with_capacity(num_workers);

        for worker_id in 0..num_workers {
            let graph = Arc::clone(&self.graph);
            let cancel = Arc::clone(&self.cancel);
            let error_tx = self.error_tx.clone();

            // Create worker state and wrap it so the main thread can retrieve stats.
            let ws = Arc::new(Mutex::new(Some(WorkerState::new(worker_id, num_workers))));
            let ws_clone = Arc::clone(&ws);
            worker_states.push(ws);

            handles.push(thread::spawn(move || {
                // Take the WorkerState out — the worker owns it for its lifetime.
                let mut state = ws_clone.lock().expect("worker state lock").take().unwrap();

                // Outer catch_unwind so panics in any stage are caught.
                let result = std::panic::catch_unwind(AssertUnwindSafe(|| {
                    worker_loop(&graph, &mut state, &cancel, &error_tx);
                }));

                if let Err(payload) = result {
                    let msg = panic_message(&payload);
                    let _ = error_tx.send(anyhow!("worker thread panicked: {msg}"));
                    cancel.store(true, Ordering::Release);
                }

                // Put the state back for stats collection.
                *ws_clone.lock().expect("worker state lock") = Some(state);
            }));
        }

        (handles, worker_states)
    }

    /// Check for errors. Returns the first error if any worker has failed.
    #[must_use]
    pub fn check_error(&self) -> Option<anyhow::Error> {
        self.error_rx.try_recv().ok()
    }

    /// Signal all workers to stop.
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Release);
    }

    /// Returns `true` if the scheduler has been cancelled.
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.cancel.load(Ordering::Acquire)
    }

    /// Return a clone of the cancellation flag.
    #[must_use]
    pub fn cancel_flag(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.cancel)
    }

    /// Return a reference to the stage graph for the calling thread to participate.
    #[must_use]
    pub fn graph(&self) -> &Arc<StageGraph> {
        &self.graph
    }

    /// Wait for all worker threads to finish, log stats, and return the first error if any.
    ///
    /// # Errors
    ///
    /// Returns an error if any worker thread panicked or sent an error.
    ///
    /// # Panics
    ///
    /// Panics if a worker-state mutex is poisoned.
    pub fn wait(
        self,
        handles: Vec<thread::JoinHandle<()>>,
        worker_states: &[WorkerStateSlot],
    ) -> Result<()> {
        for handle in handles {
            if let Err(panic_payload) = handle.join() {
                let msg = panic_message(&panic_payload);
                return Err(anyhow!("worker thread panicked: {msg}"));
            }
        }

        // Log per-worker stats.
        let num_workers = worker_states.len();
        log::info!("Pipeline worker stats ({num_workers} workers):");
        for ws in worker_states {
            if let Some(state) = ws.lock().expect("worker state lock").as_ref() {
                state.log_stats();
            }
        }

        // Drain all errors; return the first, log the rest.
        let mut first_error = None;
        while let Ok(err) = self.error_rx.try_recv() {
            if first_error.is_none() {
                first_error = Some(err);
            } else {
                log::warn!("Additional pipeline error (suppressed): {err:#}");
            }
        }
        if let Some(err) = first_error {
            return Err(err);
        }
        Ok(())
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "test indices are small enough that truncation/wrap/sign-loss cannot occur"
)]
mod tests {
    use std::sync::atomic::AtomicUsize;
    use std::time::Duration;

    use super::*;

    /// A mock stage that doubles its input.
    struct DoubleStage;

    impl PipelineStage for DoubleStage {
        type Input = i32;
        type Output = i32;

        fn process(&self, input: i32) -> Result<i32> {
            Ok(input * 2)
        }

        fn output_memory_estimate(&self, _output: &i32) -> usize {
            std::mem::size_of::<i32>()
        }
    }

    /// A mock stage that converts i32 to String.
    struct ToStringStage;

    impl PipelineStage for ToStringStage {
        type Input = i32;
        type Output = String;

        fn process(&self, input: i32) -> Result<String> {
            Ok(input.to_string())
        }

        fn output_memory_estimate(&self, output: &String) -> usize {
            output.len()
        }
    }

    /// A mock stage that returns an error on a specific input value.
    struct FailOnStage {
        fail_value: i32,
    }

    impl PipelineStage for FailOnStage {
        type Input = i32;
        type Output = i32;

        fn process(&self, input: i32) -> Result<i32> {
            if input == self.fail_value {
                Err(anyhow!("intentional failure on value {input}"))
            } else {
                Ok(input)
            }
        }

        fn output_memory_estimate(&self, _output: &i32) -> usize {
            std::mem::size_of::<i32>()
        }
    }

    /// A fast identity stage for testing priority ordering.
    struct IdentityStage;

    impl PipelineStage for IdentityStage {
        type Input = i32;
        type Output = i32;

        fn process(&self, input: i32) -> Result<i32> {
            Ok(input)
        }

        fn output_memory_estimate(&self, _output: &i32) -> usize {
            std::mem::size_of::<i32>()
        }
    }

    /// A slow stage that sleeps briefly to simulate expensive work.
    struct SlowStage {
        delay: Duration,
        call_count: AtomicUsize,
    }

    impl SlowStage {
        fn new(delay: Duration) -> Self {
            Self { delay, call_count: AtomicUsize::new(0) }
        }
    }

    impl PipelineStage for SlowStage {
        type Input = i32;
        type Output = i32;

        fn process(&self, input: i32) -> Result<i32> {
            self.call_count.fetch_add(1, Ordering::Relaxed);
            thread::sleep(self.delay);
            Ok(input)
        }

        fn output_memory_estimate(&self, _output: &i32) -> usize {
            std::mem::size_of::<i32>()
        }
    }

    /// Helper: spawn 2 worker threads running a drain-first i32->i32->O two-stage pipeline.
    fn spawn_two_stage_workers<O: Send + 'static>(
        stage1: &Arc<dyn PipelineStage<Input = i32, Output = i32>>,
        stage2: &Arc<dyn PipelineStage<Input = i32, Output = O>>,
        q_input: &WorkQueue<SequencedItem<i32>>,
        q_middle: &WorkQueue<SequencedItem<i32>>,
        q_output: &WorkQueue<SequencedItem<O>>,
        cancel: &Arc<AtomicBool>,
    ) -> Vec<thread::JoinHandle<()>> {
        let (error_tx, _error_rx) = crossbeam_channel::unbounded::<anyhow::Error>();

        (0..2)
            .map(|_| {
                let s1 = Arc::clone(stage1);
                let s2 = Arc::clone(stage2);
                let qi = q_input.clone();
                let qm = q_middle.clone();
                let qo = q_output.clone();
                let c = Arc::clone(cancel);
                let etx = error_tx.clone();

                thread::spawn(move || {
                    loop {
                        if c.load(Ordering::Acquire) {
                            break;
                        }
                        let mut did_work = false;

                        // Drain-first: try downstream (s2) first.
                        if let Some(input) = qm.try_pop() {
                            let input_mem = input.memory_estimate;
                            match s2.process(input.item) {
                                Ok(result) => {
                                    let output_mem = s2.output_memory_estimate(&result);
                                    let out = SequencedItem::new(input.seq, result, output_mem);
                                    qo.push_blocking(out, output_mem).unwrap();
                                    qm.complete_in_flight();
                                    qm.release_memory(input_mem);
                                    did_work = true;
                                }
                                Err(e) => {
                                    qm.complete_in_flight();
                                    let _ = etx.send(e);
                                    c.store(true, Ordering::Release);
                                    return;
                                }
                            }
                        } else if qm.is_closed_and_empty() {
                            qo.close();
                        }

                        if !did_work {
                            if let Some(input) = qi.try_pop() {
                                let input_mem = input.memory_estimate;
                                match s1.process(input.item) {
                                    Ok(result) => {
                                        let output_mem = s1.output_memory_estimate(&result);
                                        let out = SequencedItem::new(input.seq, result, output_mem);
                                        qm.push_blocking(out, output_mem).unwrap();
                                        qi.complete_in_flight();
                                        qi.release_memory(input_mem);
                                        did_work = true;
                                    }
                                    Err(e) => {
                                        qi.complete_in_flight();
                                        let _ = etx.send(e);
                                        c.store(true, Ordering::Release);
                                        return;
                                    }
                                }
                            } else if qi.is_closed_and_empty() {
                                qm.close();
                            }
                        }

                        if !did_work {
                            thread::sleep(Duration::from_millis(1));
                        }
                    }
                })
            })
            .collect()
    }

    #[test]
    fn test_scheduler_basic() {
        let q_input: WorkQueue<SequencedItem<i32>> = WorkQueue::new(64, usize::MAX, "input");
        let q_middle: WorkQueue<SequencedItem<i32>> = WorkQueue::new(64, usize::MAX, "middle");
        let q_output: WorkQueue<SequencedItem<String>> = WorkQueue::new(64, usize::MAX, "output");

        let cancel = Arc::new(AtomicBool::new(false));

        let num_items = 20_usize;
        for i in 0..num_items {
            let val = i as i32;
            q_input
                .push_blocking(
                    SequencedItem::new(i as u64, val, std::mem::size_of::<i32>()),
                    std::mem::size_of::<i32>(),
                )
                .unwrap();
        }
        q_input.close();

        let stage1: Arc<dyn PipelineStage<Input = i32, Output = i32>> = Arc::new(DoubleStage);
        let stage2: Arc<dyn PipelineStage<Input = i32, Output = String>> = Arc::new(ToStringStage);
        let handles =
            spawn_two_stage_workers(&stage1, &stage2, &q_input, &q_middle, &q_output, &cancel);

        let mut results: Vec<(u64, String)> = Vec::new();
        let deadline = std::time::Instant::now() + Duration::from_secs(5);
        while results.len() < num_items {
            assert!(std::time::Instant::now() < deadline, "timed out waiting for results");
            if let Some(item) = q_output.try_pop() {
                results.push((item.seq, item.item));
            } else {
                thread::sleep(Duration::from_millis(1));
            }
        }

        cancel.store(true, Ordering::Release);
        for handle in handles {
            handle.join().unwrap();
        }

        results.sort_by_key(|(seq, _)| *seq);
        for (i, (_seq, value)) in results.iter().enumerate() {
            let expected = (i as i32 * 2).to_string();
            assert_eq!(value, &expected, "seq={i}");
        }
    }

    #[test]
    fn test_scheduler_error_propagation() {
        let q_input: WorkQueue<SequencedItem<i32>> = WorkQueue::new(64, usize::MAX, "input");
        let q_output: WorkQueue<SequencedItem<i32>> = WorkQueue::new(64, usize::MAX, "output");

        let cancel = Arc::new(AtomicBool::new(false));

        for i in 0..10_i32 {
            q_input
                .push_blocking(
                    SequencedItem::new(i as u64, i, std::mem::size_of::<i32>()),
                    std::mem::size_of::<i32>(),
                )
                .unwrap();
        }
        q_input.close();

        let fail_stage = Arc::new(FailOnStage { fail_value: 5 });

        let (error_tx, error_rx) = crossbeam_channel::unbounded();
        let c = Arc::clone(&cancel);
        let qi = q_input.clone();
        let qo = q_output.clone();
        let stage = Arc::clone(&fail_stage);
        let etx = error_tx.clone();

        let handle = thread::spawn(move || {
            loop {
                if c.load(Ordering::Acquire) {
                    break;
                }
                if let Some(input) = qi.try_pop() {
                    let input_mem = input.memory_estimate;
                    match stage.process(input.item) {
                        Ok(result) => {
                            let output_mem = stage.output_memory_estimate(&result);
                            let out = SequencedItem::new(input.seq, result, output_mem);
                            qo.push_blocking(out, output_mem).unwrap();
                            qi.complete_in_flight();
                            qi.release_memory(input_mem);
                        }
                        Err(e) => {
                            qi.complete_in_flight();
                            let _ = etx.send(e);
                            c.store(true, Ordering::Release);
                            return;
                        }
                    }
                } else if qi.is_closed_and_empty() {
                    break;
                } else {
                    thread::sleep(Duration::from_millis(1));
                }
            }
        });

        handle.join().unwrap();

        let err = error_rx.try_recv().expect("should have received an error");
        let err_msg = err.to_string();
        assert!(err_msg.contains("intentional failure"), "unexpected error message: {err_msg}");
    }

    #[test]
    fn test_scheduler_drain_first() {
        let q_input: WorkQueue<SequencedItem<i32>> = WorkQueue::new(128, usize::MAX, "input");
        let q_middle: WorkQueue<SequencedItem<i32>> = WorkQueue::new(128, usize::MAX, "middle");
        let q_output: WorkQueue<SequencedItem<i32>> = WorkQueue::new(128, usize::MAX, "output");

        let cancel = Arc::new(AtomicBool::new(false));

        let slow = Arc::new(SlowStage::new(Duration::from_millis(5)));
        let fast = Arc::new(IdentityStage);

        let num_items = 40_usize;
        for i in 0..num_items {
            q_input
                .push_blocking(
                    SequencedItem::new(i as u64, i as i32, std::mem::size_of::<i32>()),
                    std::mem::size_of::<i32>(),
                )
                .unwrap();
        }
        q_input.close();

        let fast_stage: Arc<dyn PipelineStage<Input = i32, Output = i32>> = fast;
        let slow_stage: Arc<dyn PipelineStage<Input = i32, Output = i32>> = slow;
        let handles = spawn_two_stage_workers(
            &fast_stage,
            &slow_stage,
            &q_input,
            &q_middle,
            &q_output,
            &cancel,
        );

        let mut max_middle_len: usize = 0;
        let mut collected: usize = 0;
        let deadline = std::time::Instant::now() + Duration::from_secs(10);

        while collected < num_items {
            assert!(std::time::Instant::now() < deadline, "timed out waiting for results");
            let middle_len = q_middle.len();
            if middle_len > max_middle_len {
                max_middle_len = middle_len;
            }
            if q_output.try_pop().is_some() {
                collected += 1;
            } else {
                thread::sleep(Duration::from_millis(1));
            }
        }

        cancel.store(true, Ordering::Release);
        for handle in handles {
            handle.join().unwrap();
        }

        assert!(
            max_middle_len < num_items,
            "middle queue grew to {max_middle_len} (total items: {num_items}); \
             drain-first priority should prevent this"
        );
    }

    #[test]
    fn test_write_preferred_starts_on_compress() {
        // WritePreferred starts with current_step = Compress.
        let mut sched = PipelineScheduler::new(3, 4);
        assert_eq!(sched.role, WorkerRole::WritePreferred);
        let bp = BackpressureState { output_high: false, memory_high: false };
        let priorities = sched.get_priorities(&bp);
        // Write is the exclusive role step, so it goes first.
        assert_eq!(priorities[0], PipelineStep::Write);
        // Compress is the current_step (sticky start), so it goes second.
        assert_eq!(priorities[1], PipelineStep::Compress);
    }

    #[test]
    fn test_consensus_worker_starts_on_consensus() {
        let mut sched = PipelineScheduler::new(0, 4);
        assert_eq!(sched.role, WorkerRole::ConsensusWorker);
        let bp = BackpressureState { output_high: false, memory_high: false };
        let priorities = sched.get_priorities(&bp);
        // ConsensusWorker has no exclusive role, current_step = Consensus.
        assert_eq!(priorities[0], PipelineStep::Consensus);
    }

    #[test]
    fn test_sticky_on_success() {
        let mut sched = PipelineScheduler::new(0, 4);
        // Record success on Filter — current_step should become Filter.
        sched.record_outcome(PipelineStep::Filter, true);
        let bp = BackpressureState { output_high: false, memory_high: false };
        let priorities = sched.get_priorities(&bp);
        // Filter is now the current step (sticky), so it goes first.
        assert_eq!(priorities[0], PipelineStep::Filter);
    }

    #[test]
    fn test_write_preferred_pivots_after_write() {
        let mut sched = PipelineScheduler::new(3, 4);
        assert_eq!(sched.role, WorkerRole::WritePreferred);
        // After successful Write, WritePreferred pivots to Compress.
        sched.record_outcome(PipelineStep::Write, true);
        let bp = BackpressureState { output_high: false, memory_high: false };
        let priorities = sched.get_priorities(&bp);
        assert_eq!(priorities[0], PipelineStep::Write); // exclusive role
        assert_eq!(priorities[1], PipelineStep::Compress); // pivoted current_step
    }

    #[test]
    fn test_group_preferred_pivots_after_group() {
        let mut sched = PipelineScheduler::new(1, 4);
        assert_eq!(sched.role, WorkerRole::GroupAssignPreferred);
        // After successful GroupAssign, GroupAssignPreferred pivots to Consensus.
        sched.record_outcome(PipelineStep::GroupAssign, true);
        let bp = BackpressureState { output_high: false, memory_high: false };
        let priorities = sched.get_priorities(&bp);
        assert_eq!(priorities[0], PipelineStep::GroupAssign); // exclusive role
        assert_eq!(priorities[1], PipelineStep::Consensus); // pivoted current_step
    }

    #[test]
    fn test_backpressure_moves_compress_up() {
        let mut sched = PipelineScheduler::new(0, 4);
        assert_eq!(sched.role, WorkerRole::ConsensusWorker);
        let bp = BackpressureState { output_high: true, memory_high: false };
        let priorities = sched.get_priorities(&bp);
        // With output_high, Compress should be early (before GroupAssign).
        let compress_pos = priorities.iter().position(|s| *s == PipelineStep::Compress).unwrap();
        let group_pos = priorities.iter().position(|s| *s == PipelineStep::GroupAssign).unwrap();
        assert!(
            compress_pos < group_pos,
            "Compress should have higher priority than GroupAssign under backpressure"
        );
    }

    #[test]
    fn test_role_assignment() {
        // 4 workers: 0=ConsensusWorker, 1=GroupAssignPreferred, 2=ConsensusWorker, 3=WritePreferred
        assert_eq!(PipelineScheduler::new(0, 4).role, WorkerRole::ConsensusWorker);
        assert_eq!(PipelineScheduler::new(1, 4).role, WorkerRole::GroupAssignPreferred);
        assert_eq!(PipelineScheduler::new(2, 4).role, WorkerRole::ConsensusWorker);
        assert_eq!(PipelineScheduler::new(3, 4).role, WorkerRole::WritePreferred);

        // 3 workers: 0=ConsensusWorker, 1=ConsensusWorker (not enough for GroupAssign), 2=WritePreferred
        assert_eq!(PipelineScheduler::new(0, 3).role, WorkerRole::ConsensusWorker);
        assert_eq!(PipelineScheduler::new(1, 3).role, WorkerRole::ConsensusWorker);
        assert_eq!(PipelineScheduler::new(2, 3).role, WorkerRole::WritePreferred);

        // 2 workers: 0=ConsensusWorker, 1=WritePreferred
        assert_eq!(PipelineScheduler::new(0, 2).role, WorkerRole::ConsensusWorker);
        assert_eq!(PipelineScheduler::new(1, 2).role, WorkerRole::WritePreferred);

        // 1 worker: 0=WritePreferred (it's the last worker)
        assert_eq!(PipelineScheduler::new(0, 1).role, WorkerRole::WritePreferred);
    }

    #[test]
    fn test_worker_state_held_items() {
        let mut ws = WorkerState::new(0, 1);
        assert!(!ws.has_any_held_items());

        ws.held_compress = Some((SequencedItem::new(0, CompressedBatch { blocks: vec![] }, 0), 0));
        assert!(ws.has_any_held_items());
    }

    #[test]
    fn test_pipeline_step_from_index() {
        assert_eq!(PipelineStep::from_index(0), PipelineStep::GroupAssign);
        assert_eq!(PipelineStep::from_index(1), PipelineStep::Consensus);
        assert_eq!(PipelineStep::from_index(2), PipelineStep::Filter);
        assert_eq!(PipelineStep::from_index(3), PipelineStep::Compress);
        assert_eq!(PipelineStep::from_index(4), PipelineStep::Write);
    }
}

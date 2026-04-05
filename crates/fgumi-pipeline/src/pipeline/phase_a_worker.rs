//! Phase A worker state and pool execution loop.
//!
//! Each pool worker owns a [`PhaseAWorkerState`] with held-item slots and scheduler state.
//! The [`start_phase_a_pool`] function spawns worker threads that execute the Phase A loop
//! until all queues are drained or cancellation is signaled.

use std::panic::AssertUnwindSafe;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

use anyhow::anyhow;
use crossbeam_channel::{Receiver, Sender};
use lru::LruCache;

use fgumi_lib::correct::UmiMatch;

use crate::pipeline::phase_a::{PhaseAScheduler, PhaseAStep};
use crate::pipeline::phase_a_graph::{OrdinalPair, PhaseAGraph};
use crate::pipeline::phase_a_steps;
use crate::pipeline::scheduler::StepResult;
use crate::pipeline::worker_util::{WorkerBackoff, panic_message};

// ============================================================================
// Worker stats
// ============================================================================

/// Per-worker statistics for Phase A diagnostics.
pub struct PhaseAWorkerStats {
    /// Per-stage attempt count.
    pub stage_attempts: [AtomicU64; PhaseAStep::NUM_STEPS],
    /// Per-stage success count.
    pub stage_successes: [AtomicU64; PhaseAStep::NUM_STEPS],
    /// Number of idle cycles (no work found in any stage).
    pub idle_cycles: AtomicU64,
    /// Cumulative nanoseconds spent in backoff sleep.
    pub idle_nanos: AtomicU64,
}

impl PhaseAWorkerStats {
    /// Create a new zeroed stats instance.
    fn new() -> Self {
        Self {
            stage_attempts: std::array::from_fn(|_| AtomicU64::new(0)),
            stage_successes: std::array::from_fn(|_| AtomicU64::new(0)),
            idle_cycles: AtomicU64::new(0),
            idle_nanos: AtomicU64::new(0),
        }
    }
}

// ============================================================================
// PhaseAWorkerState
// ============================================================================

/// Per-worker mutable state for Phase A pool stages.
///
/// Each worker thread owns one instance. Held items are stored here when output queues
/// are full, enabling non-blocking push semantics.
pub struct PhaseAWorkerState {
    /// Worker ID (0-indexed).
    pub worker_id: usize,

    // ---- Extract held item ----
    /// Held ordinal-tagged `(ordinal, Template, fastq_bytes)` triple from Extract
    /// (`q_extracted` was full).
    pub held_extract_pair: Option<OrdinalPair>,
    /// Memory estimate for the held extract pair.
    pub held_extract_pair_mem: usize,

    // ---- Correct held item ----
    /// Held corrected ordinal-tagged `(ordinal, Template, fastq_bytes)` triple
    /// (`q_corrected` was full).
    pub held_correct_pair: Option<OrdinalPair>,
    /// Memory estimate for the held correct pair.
    pub held_correct_pair_mem: usize,
    /// Per-worker LRU cache for UMI segment matching.
    /// `None` when correction is disabled or `cache_size` is 0.
    pub umi_cache: Option<LruCache<Vec<u8>, UmiMatch>>,

    // ---- Scheduler and backoff ----
    /// Phase A scheduler for priority determination.
    pub scheduler: PhaseAScheduler,
    /// Per-worker statistics.
    pub stats: PhaseAWorkerStats,
    /// Adaptive backoff with jitter.
    backoff: WorkerBackoff,
}

impl PhaseAWorkerState {
    /// Create a new worker state.
    ///
    /// # Panics
    ///
    /// Panics if `umi_cache_size` is non-zero but cannot be converted to a `NonZero<usize>`.
    #[must_use]
    pub fn new(
        worker_id: usize,
        active_steps: [bool; PhaseAStep::NUM_STEPS],
        umi_cache_size: usize,
    ) -> Self {
        let umi_cache = if umi_cache_size > 0 {
            Some(LruCache::new(std::num::NonZero::new(umi_cache_size).unwrap()))
        } else {
            None
        };
        Self {
            worker_id,
            held_extract_pair: None,
            held_extract_pair_mem: 0,
            held_correct_pair: None,
            held_correct_pair_mem: 0,
            umi_cache,
            scheduler: PhaseAScheduler::new(active_steps),
            stats: PhaseAWorkerStats::new(),
            backoff: WorkerBackoff::new(worker_id),
        }
    }

    /// Returns `true` if this worker is holding any items.
    pub fn has_any_held_items(&self) -> bool {
        self.held_extract_pair.is_some() || self.held_correct_pair.is_some()
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

    /// Log stats for this worker.
    #[expect(clippy::cast_precision_loss, reason = "nanos fits in f64 for display")]
    pub fn log_stats(&self) {
        log::info!("  Phase A Worker {}:", self.worker_id);
        for &step in &PhaseAStep::ALL {
            let idx = step as usize;
            let attempts = self.stats.stage_attempts[idx].load(Ordering::Relaxed);
            let successes = self.stats.stage_successes[idx].load(Ordering::Relaxed);
            if attempts > 0 {
                log::info!(
                    "    {:<14} {attempts:>8} attempts, {successes:>8} success",
                    step.name()
                );
            }
        }
        let idle_cycles = self.stats.idle_cycles.load(Ordering::Relaxed);
        let idle_nanos = self.stats.idle_nanos.load(Ordering::Relaxed);
        let idle_secs = idle_nanos as f64 / 1_000_000_000.0;
        log::info!("    idle: {idle_cycles} cycles, {idle_secs:.3}s");
    }
}

// ============================================================================
// Worker loop
// ============================================================================

/// The main loop executed by each Phase A worker thread.
///
/// Each iteration: (1) sample backpressure and get priorities, (2) try one step in
/// priority order, (3) backoff if idle.
fn phase_a_worker_loop(
    graph: &PhaseAGraph,
    worker: &mut PhaseAWorkerState,
    cancel: &AtomicBool,
    error_tx: &Sender<anyhow::Error>,
) {
    loop {
        // Exit when cancelled AND no held items remain.
        if cancel.load(Ordering::Acquire) && !worker.has_any_held_items() {
            break;
        }

        // Check if all pool work is done.
        if graph.is_pool_done() && !worker.has_any_held_items() {
            break;
        }

        let mut did_work = false;

        // 1. Sample backpressure and get priority order.
        let bp = graph.sample_backpressure();
        let priorities = worker.scheduler.get_priorities(&bp);
        // Copy to avoid borrow conflict with worker mutation in the loop.
        let priorities_copy: Vec<PhaseAStep> = priorities.to_vec();

        // 2. Execute one stage in priority order.
        for &step in &priorities_copy {
            if cancel.load(Ordering::Acquire) {
                break;
            }

            let idx = step as usize;
            worker.stats.stage_attempts[idx].fetch_add(1, Ordering::Relaxed);

            let result = match step {
                PhaseAStep::Extract => phase_a_steps::try_step_extract(graph, worker),
                PhaseAStep::Correct => phase_a_steps::try_step_correct(graph, worker),
                PhaseAStep::ToFastq => phase_a_steps::try_step_to_fastq(graph, worker),
            };

            match result {
                StepResult::Success => {
                    worker.stats.stage_successes[idx].fetch_add(1, Ordering::Relaxed);
                    worker.scheduler.record_outcome(step, true);
                    did_work = true;
                    break;
                }
                StepResult::OutputFull => {
                    worker.scheduler.record_outcome(step, false);
                }
                StepResult::InputEmpty => {}
                StepResult::Error(e) => {
                    let _ = error_tx.send(e);
                    cancel.store(true, Ordering::Release);
                    return;
                }
            }
        }

        // 3. Adaptive backoff if no work done.
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
// Pool execution
// ============================================================================

/// Shared slot for returning worker state after thread completion.
pub type PhaseAWorkerSlot = Arc<Mutex<Option<PhaseAWorkerState>>>;

/// Run the Phase A worker pool.
///
/// Spawns `num_workers` threads that execute the Phase A work-stealing loop. Returns
/// the join handles and worker state slots for stats collection.
///
/// The active steps determine which stages are enabled. For the initial implementation
/// (no Correct), the active steps are: Extract, `ToFastq`.
///
/// # Panics
///
/// Panics if a worker-state mutex is poisoned or thread spawning fails.
#[must_use]
pub fn start_phase_a_pool(
    graph: &Arc<PhaseAGraph>,
    num_workers: usize,
    active_steps: [bool; PhaseAStep::NUM_STEPS],
    cancel: &Arc<AtomicBool>,
    umi_cache_size: usize,
) -> (Vec<thread::JoinHandle<()>>, Vec<PhaseAWorkerSlot>, Receiver<anyhow::Error>) {
    let (error_tx, error_rx) = crossbeam_channel::unbounded();
    let mut handles = Vec::with_capacity(num_workers);
    let mut worker_slots = Vec::with_capacity(num_workers);

    for worker_id in 0..num_workers {
        let graph = Arc::clone(graph);
        let cancel = Arc::clone(cancel);
        let error_tx = error_tx.clone();

        let ws = Arc::new(Mutex::new(Some(PhaseAWorkerState::new(
            worker_id,
            active_steps,
            umi_cache_size,
        ))));
        let ws_clone = Arc::clone(&ws);
        worker_slots.push(ws);

        handles.push(
            thread::Builder::new()
                .name(format!("phase-a-{worker_id}"))
                .spawn(move || {
                    let mut state = ws_clone.lock().expect("worker state lock").take().unwrap();

                    let result = std::panic::catch_unwind(AssertUnwindSafe(|| {
                        phase_a_worker_loop(&graph, &mut state, &cancel, &error_tx);
                    }));

                    if let Err(payload) = result {
                        let msg = panic_message(&payload);
                        let _ = error_tx.send(anyhow!("Phase A worker thread panicked: {msg}"));
                        cancel.store(true, Ordering::Release);
                    }

                    // Put the state back for stats collection.
                    *ws_clone.lock().expect("worker state lock") = Some(state);
                })
                .expect("failed to spawn Phase A worker thread"),
        );
    }

    (handles, worker_slots, error_rx)
}

/// Wait for all Phase A worker threads to finish and return the first error if any.
///
/// # Errors
///
/// Returns an error if any worker thread panicked or sent an error.
///
/// # Panics
///
/// Panics if a worker-state mutex is poisoned.
pub fn wait_phase_a_pool(
    handles: Vec<thread::JoinHandle<()>>,
    worker_slots: &[PhaseAWorkerSlot],
    error_rx: &Receiver<anyhow::Error>,
) -> anyhow::Result<()> {
    for handle in handles {
        if let Err(payload) = handle.join() {
            let msg = panic_message(&payload);
            return Err(anyhow!("Phase A worker thread panicked: {msg}"));
        }
    }

    // Log worker stats.
    log::info!("Phase A worker stats:");
    for slot in worker_slots {
        if let Some(ref state) = *slot.lock().expect("worker state lock") {
            state.log_stats();
        }
    }

    // Check for errors.
    if let Ok(e) = error_rx.try_recv() {
        return Err(e);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_worker_state_initial() {
        let active = [true; PhaseAStep::NUM_STEPS];
        let worker = PhaseAWorkerState::new(0, active, 0);
        assert_eq!(worker.worker_id, 0);
        assert!(!worker.has_any_held_items());
    }

    #[test]
    fn test_worker_state_held_items() {
        use fgumi_lib::template::Template;
        let active = [true; PhaseAStep::NUM_STEPS];
        let mut worker = PhaseAWorkerState::new(0, active, 0);
        assert!(!worker.has_any_held_items());

        worker.held_extract_pair = Some((0, Template::new(b"test".to_vec()), vec![1, 2, 3]));
        assert!(worker.has_any_held_items());
    }

    #[test]
    fn test_worker_state_with_cache() {
        let active = [true; PhaseAStep::NUM_STEPS];
        let worker = PhaseAWorkerState::new(0, active, 100);
        assert!(worker.umi_cache.is_some());
        assert!(!worker.has_any_held_items());
    }

    #[test]
    fn test_worker_state_no_cache() {
        let active = [true; PhaseAStep::NUM_STEPS];
        let worker = PhaseAWorkerState::new(0, active, 0);
        assert!(worker.umi_cache.is_none());
    }

    #[test]
    fn test_worker_state_held_correct_pair() {
        use fgumi_lib::template::Template;
        let active = [true; PhaseAStep::NUM_STEPS];
        let mut worker = PhaseAWorkerState::new(0, active, 0);
        assert!(!worker.has_any_held_items());

        worker.held_correct_pair = Some((0, Template::new(b"test".to_vec()), vec![1, 2, 3]));
        assert!(worker.has_any_held_items());
    }
}

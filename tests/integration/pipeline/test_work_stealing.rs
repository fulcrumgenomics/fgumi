//! Gate tests that make it structurally impossible to regress to
//! one-thread-per-stage. See docs/superpowers/plans/2026-04-13-work-stealing-pool.md.
//!
//! These tests require the work-stealing pool implementation; they will FAIL
//! on any per-stage-pinned worker model.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use fgumi_lib::runall::engine::cancel::CancelToken;
use fgumi_lib::runall::engine::driver::{PipelineConfig, StageEntry, run_multi_stage};
use fgumi_lib::runall::engine::stage::{Parallelism, Stage};
use fgumi_lib::runall::engine::stage_source::{ErasedStageSource, StageSource};

use super::mock_stages::{AddStage, CollectSink, RangeSource};

use crate::helpers::thread_registry;

/// Assert that every live thread name is a well-known pipeline role and that
/// exactly `worker_threads` of them are pool workers.
fn assert_only_pool_and_expected_threads_exist(worker_threads: usize) {
    let registry = thread_registry::snapshot();
    let pool_workers: Vec<_> = registry.iter().filter(|n| n.starts_with("pool_worker_")).collect();
    assert_eq!(
        pool_workers.len(),
        worker_threads,
        "expected exactly {worker_threads} pool workers, got {pool_workers:?} \
         (full registry: {registry:?})",
    );
    for name in &registry {
        assert!(
            name.starts_with("pool_worker_")
                || name == "pipeline_source"
                || name == "pipeline_sink"
                || name == "pipeline_watchdog",
            "unexpected thread name {name:?}; plan 2b-1 forbids per-stage worker threads",
        );
    }
}

#[test]
fn test_a_pool_has_exactly_n_workers() {
    thread_registry::clear();

    let out = Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out });
    let source = Box::new(RangeSource { n: 100 });

    let stages: Vec<StageEntry> = vec![
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
    ];

    let worker_threads = 4;
    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads,
            queue_capacity: 64,
            queue_memory_limit: 1024 * 1024,
            global_memory_limit: 16 * 1024 * 1024,
        },
    )
    .unwrap();

    assert_only_pool_and_expected_threads_exist(worker_threads);
}

// ----------------------------------------------------------------------------
// Test B: any worker can run any stage (not pinned to a single stage).
// ----------------------------------------------------------------------------

/// A stage that records which worker-thread-id executed it.
#[derive(Clone)]
struct TrackedStage {
    stage_id: u32,
    log: Arc<std::sync::Mutex<Vec<(std::thread::ThreadId, u32)>>>,
}

impl Stage for TrackedStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        self.log.lock().unwrap().push((std::thread::current().id(), self.stage_id));
        out(input + u64::from(self.stage_id));
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, _o: &u64) -> usize {
        8
    }

    fn name(&self) -> &'static str {
        "TrackedStage"
    }
}

#[test]
fn test_b_any_worker_runs_any_stage() {
    // Skip under coverage instrumentation: per-function counters slow hot
    // loops enough that the scheduler hands out items faster than workers
    // can steal across stages, and one worker ends up pinned to a single
    // stage in the log. Not a real regression of the work-stealing scheduler.
    if std::env::var_os("CARGO_LLVM_COV").is_some() {
        eprintln!("SKIP: work-stealing distribution test is timing-sensitive under coverage");
        return;
    }

    let log: Arc<std::sync::Mutex<Vec<(std::thread::ThreadId, u32)>>> =
        Arc::new(std::sync::Mutex::new(Vec::new()));

    let out = Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out });
    let source = Box::new(RangeSource { n: 500 });

    let stages: Vec<StageEntry> = (0..3u32)
        .map(|i| {
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(TrackedStage {
                stage_id: i,
                log: log.clone(),
            })))
        })
        .collect();

    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads: 4,
            queue_capacity: 32,
            queue_memory_limit: 1024 * 1024,
            global_memory_limit: 16 * 1024 * 1024,
        },
    )
    .unwrap();

    // Build: worker_tid -> set of stage_ids it executed.
    let log = log.lock().unwrap();
    let mut worker_stages: std::collections::HashMap<
        std::thread::ThreadId,
        std::collections::HashSet<u32>,
    > = std::collections::HashMap::new();
    for (tid, sid) in log.iter() {
        worker_stages.entry(*tid).or_default().insert(*sid);
    }

    // Must see at least 2 distinct workers (with 4 workers and 500 items,
    // at least 2 are guaranteed to participate).
    assert!(
        worker_stages.len() >= 2,
        "expected >=2 workers to participate, saw {}",
        worker_stages.len()
    );

    // Every participating worker must have run at least 2 distinct stages.
    // With one-thread-per-stage, each worker runs exactly 1 stage — this fails.
    for (tid, stages_for_worker) in &worker_stages {
        assert!(
            stages_for_worker.len() >= 2,
            "worker {:?} only ran {} stages: {:?} — pinning to one stage detected",
            tid,
            stages_for_worker.len(),
            stages_for_worker,
        );
    }
}

// ----------------------------------------------------------------------------
// Test C: workers cluster on a 10x-slower stage.
// ----------------------------------------------------------------------------

#[derive(Clone)]
struct SlowStage {
    slow: bool,
    /// Per-stage in-flight counter; used to track how many workers are
    /// concurrently inside `process()` for this stage.
    in_flight: Arc<AtomicUsize>,
    /// Per-stage max-observed-concurrency; updated under `in_flight`.
    max_concurrent: Arc<AtomicUsize>,
}

impl Stage for SlowStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        let current = self.in_flight.fetch_add(1, Ordering::SeqCst) + 1;
        self.max_concurrent.fetch_max(current, Ordering::SeqCst);
        if self.slow {
            // Simulate 10x slower work.
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        self.in_flight.fetch_sub(1, Ordering::SeqCst);
        out(input);
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, _o: &u64) -> usize {
        8
    }

    fn name(&self) -> &'static str {
        "SlowStage"
    }
}

#[test]
fn test_c_workers_cluster_on_slow_stage() {
    // Skip under coverage instrumentation: per-function counter overhead
    // distorts the timing that the "workers cluster on slow stage" signal
    // relies on. Same rationale as test_b / test_e.
    if std::env::var_os("CARGO_LLVM_COV").is_some() {
        eprintln!("SKIP: worker-clustering test is timing-sensitive under coverage");
        return;
    }

    // For each stage, track how many workers are in `process()` at any
    // given moment. The slow middle stage is the bottleneck; the scheduler
    // should direct MULTIPLE workers to run it concurrently (clustering)
    // while the fast stages run with at most 1 concurrent worker at a time.
    let in_flight: Vec<Arc<AtomicUsize>> = (0..3).map(|_| Arc::new(AtomicUsize::new(0))).collect();
    let max_concurrent: Vec<Arc<AtomicUsize>> =
        (0..3).map(|_| Arc::new(AtomicUsize::new(0))).collect();

    let out = Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out });
    let source = Box::new(RangeSource { n: 200 });

    // Use the work-stealing entry point so each Parallel stage gets
    // per-worker instances — otherwise the pool coerces the single instance
    // to Sequential and the slow-stage clustering signal is invisible.
    let stages: Vec<StageEntry> = (0..3usize)
        .map(|i| {
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(SlowStage {
                slow: i == 1,
                in_flight: in_flight[i].clone(),
                max_concurrent: max_concurrent[i].clone(),
            })))
        })
        .collect();

    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads: 4,
            queue_capacity: 16,
            queue_memory_limit: 1024 * 1024,
            global_memory_limit: 16 * 1024 * 1024,
        },
    )
    .unwrap();

    let fast_peak_0 = max_concurrent[0].load(Ordering::SeqCst);
    let slow_peak = max_concurrent[1].load(Ordering::SeqCst);
    let fast_peak_2 = max_concurrent[2].load(Ordering::SeqCst);

    // Clustering: the slow stage should be entered by multiple workers
    // concurrently. With 4 workers we expect at least 2 to run stage 1 at
    // once. Under per-stage pinning this would be 1 — the stage has only
    // one worker — so the test fails.
    assert!(
        slow_peak >= 2,
        "expected workers to cluster on slow stage; peak concurrency was \
         {slow_peak} (fast 0: {fast_peak_0}, fast 2: {fast_peak_2})",
    );
}

// ----------------------------------------------------------------------------
// Test D: a Sequential stage never has two workers in process() concurrently.
// ----------------------------------------------------------------------------

#[derive(Clone)]
struct SequentialCheckedStage {
    in_flight: Arc<AtomicUsize>,
    max_observed: Arc<AtomicUsize>,
}

impl Stage for SequentialCheckedStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        let current = self.in_flight.fetch_add(1, Ordering::SeqCst) + 1;
        // Record the maximum observed concurrency.
        self.max_observed.fetch_max(current, Ordering::SeqCst);
        std::thread::sleep(std::time::Duration::from_micros(100));
        self.in_flight.fetch_sub(1, Ordering::SeqCst);
        out(input);
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Sequential
    }

    fn output_memory_estimate(&self, _o: &u64) -> usize {
        8
    }

    fn name(&self) -> &'static str {
        "SequentialCheckedStage"
    }
}

#[test]
fn test_d_sequential_stage_is_mutually_exclusive() {
    let in_flight = Arc::new(AtomicUsize::new(0));
    let max_observed = Arc::new(AtomicUsize::new(0));

    let out = Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out });
    let source = Box::new(RangeSource { n: 200 });

    let stages: Vec<StageEntry> = vec![StageEntry::Pool(ErasedStageSource::from_source(
        StageSource::Clone(SequentialCheckedStage {
            in_flight: in_flight.clone(),
            max_observed: max_observed.clone(),
        }),
    ))];

    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads: 8,
            queue_capacity: 32,
            queue_memory_limit: 1024 * 1024,
            global_memory_limit: 16 * 1024 * 1024,
        },
    )
    .unwrap();

    let peak = max_observed.load(Ordering::SeqCst);
    assert_eq!(peak, 1, "Sequential stage saw {peak} concurrent workers; must be 1");
}

// ----------------------------------------------------------------------------
// Test E: CPU-bound Parallel stage shows >=3x speedup from 1 -> 8 workers.
// ----------------------------------------------------------------------------

#[derive(Clone)]
struct CpuBoundStage;

impl Stage for CpuBoundStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        // Burn real CPU per item (no sleeps).
        let mut h: u64 = input;
        for _ in 0..10_000 {
            h = h.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        }
        out(h);
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, _o: &u64) -> usize {
        8
    }

    fn name(&self) -> &'static str {
        "CpuBoundStage"
    }
}

fn run_cpu_bound(worker_threads: usize, n: u64) -> std::time::Duration {
    let out = Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out });
    let source = Box::new(RangeSource { n });
    let stages: Vec<StageEntry> =
        vec![StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(CpuBoundStage)))];
    let start = std::time::Instant::now();
    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads,
            queue_capacity: 256,
            queue_memory_limit: 16 * 1024 * 1024,
            global_memory_limit: 64 * 1024 * 1024,
        },
    )
    .unwrap();
    start.elapsed()
}

#[test]
fn test_e_throughput_scales_with_workers() {
    // Skip under coverage instrumentation: llvm-cov's per-function counters
    // serialise hot loops enough that multi-worker runs measure slower than
    // single-worker, which is not a real regression of the scheduler.
    // `CARGO_LLVM_COV` is set by `cargo llvm-cov nextest`.
    if std::env::var_os("CARGO_LLVM_COV").is_some() {
        eprintln!("SKIP: throughput scaling test is meaningless under coverage");
        return;
    }
    // The per-iteration CPU work is sized so that a serial run
    // dominates scheduler/pool overhead and a parallel run shows real
    // speedup. We run each configuration several times and take the best
    // (minimum) elapsed wall-time to reduce noise from concurrent tests
    // in a CI runner — the shortest runs reflect when the OS scheduler
    // gave us the CPUs we asked for.
    //
    // Worker count is scaled to the machine's available parallelism so the
    // test runs meaningfully on 2-core GitHub Actions runners as well as
    // big-iron developer boxes. The speedup threshold scales with the
    // available cores: we expect at least 40% utilisation (e.g. 1.6× on
    // 4 workers), well above the 1.0× a broken one-thread-per-stage
    // scheduler would produce, yet tolerant of shared-runner scheduler noise.
    let available = std::thread::available_parallelism().map(std::num::NonZero::get).unwrap_or(2);
    let workers = available.max(2);
    let n = 50_000;
    let trials = 3;

    let mut best_t1 = std::time::Duration::MAX;
    for _ in 0..trials {
        let t = run_cpu_bound(1, n);
        if t < best_t1 {
            best_t1 = t;
        }
    }

    let mut best_parallel = std::time::Duration::MAX;
    for _ in 0..trials {
        let t = run_cpu_bound(workers, n);
        if t < best_parallel {
            best_parallel = t;
        }
    }

    #[allow(clippy::cast_precision_loss)]
    let expected_min_speedup = (workers as f64) * 0.4;
    let speedup = best_t1.as_secs_f64() / best_parallel.as_secs_f64();
    assert!(
        speedup >= expected_min_speedup,
        "{workers}-worker speedup {:.2}x (best t1={:.2}s, best t{workers}={:.2}s); \
         expected >={:.2}x — workers not parallelizing",
        speedup,
        best_t1.as_secs_f64(),
        best_parallel.as_secs_f64(),
        expected_min_speedup,
    );
}

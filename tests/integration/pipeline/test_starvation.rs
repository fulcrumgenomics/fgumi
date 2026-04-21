//! Regression tests for held-item starvation (v1 issues #251, #253).
//!
//! The starvation pattern: when every worker holds a blocked output item
//! (because the downstream queue is full), workers previously tight-spun on
//! retrying ONLY the held slot with `continue`, so nobody popped new input.
//! When the downstream queue eventually drained, workers were still stuck in
//! the held-only retry path.
//!
//! The fix (see `worker.rs` and `run_erased_worker` in `pipeline.rs`): after a
//! held-push fails, fall through to the normal iteration tail and only yield
//! once per iteration. That keeps cancellation observation, held-retry, and
//! new-work-pop interleaved.
//!
//! With only one worker per stage, the starvation is not directly reachable
//! (there is no "every worker holds a blocked item" state — a single worker
//! that cannot push simply cannot pop, and the downstream stage is still
//! running). The regression test is therefore kept as a documented `#[ignore]`
//! placeholder until multi-worker-per-stage is wired up. The fix itself is
//! still applied because the held-only tight-spin is also a CPU waster and
//! the code pattern is a footgun for future multi-worker expansion.

use std::sync::Arc;
use std::sync::Mutex;

use fgumi_lib::runall::engine::cancel::CancelToken;
use fgumi_lib::runall::engine::driver::{PipelineConfig, StageEntry, run_multi_stage};
use fgumi_lib::runall::engine::stage_source::{ErasedStageSource, StageSource};

use super::mock_stages::{AddStage, CollectSink, RangeSource};

/// Tight-capacity pipeline with a slow sink. Exercises the held-item path
/// heavily: the output queue between stage and sink saturates quickly, each
/// worker's push fails, the held slot fills, and the worker must still
/// observe drain and cancellation. After the F17 fix, the pipeline completes
/// successfully rather than spinning on held-only retries.
#[test]
fn test_tight_queue_completes_without_starvation() {
    let out = Arc::new(Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out: out.clone() });
    let source = Box::new(RangeSource { n: 200 });

    let stages: Vec<StageEntry> = vec![
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
    ];

    // Deliberately tight: every queue boundary will block frequently.
    let cfg = PipelineConfig {
        worker_threads: 1,
        queue_capacity: 2,
        queue_memory_limit: 64,
        global_memory_limit: 1_000_000,
    };

    run_multi_stage::<u64, u64>(source, stages, sink, &CancelToken::new(), cfg).unwrap();

    let got = out.lock().unwrap().clone();
    let expected: Vec<_> = (0..200u64).map(|i| i + 3).collect();
    assert_eq!(got, expected);
}

/// Placeholder for the true multi-worker starvation reproducer. Single-worker
/// pipelines cannot exhibit the "every worker stuck on held" pathology, so
/// this test is ignored until multi-worker-per-stage scheduling lands.
#[test]
#[ignore = "requires multi-worker-per-stage; see v1 issues #251, #253"]
#[allow(dead_code)]
fn test_multi_worker_held_starvation_reproducer() {
    // Intentionally empty. Reintroduce when multi-worker stages are wired up.
}

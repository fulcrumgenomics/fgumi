//! End-to-end tests for the `pipeline` engine using mock stages.

use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

use fgumi_lib::runall::engine::cancel::CancelToken;
use fgumi_lib::runall::engine::driver::{PipelineConfig, StageEntry, run_multi_stage};
use fgumi_lib::runall::engine::stage_source::{ErasedStageSource, StageSource};

use super::mock_stages::{AddStage, CollectSink, CountingStage, MulStage, RangeSource};

#[test]
fn test_three_stage_pipeline_correctness() {
    // Source: 0..20. +3, *2, +1. Expected: (i + 3) * 2 + 1
    let out = Arc::new(Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out: out.clone() });
    let source = Box::new(RangeSource { n: 20 });

    let stages: Vec<StageEntry> = vec![
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 3 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(MulStage { n: 2 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
    ];

    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads: 4,
            queue_capacity: 8,
            queue_memory_limit: 4096,
            global_memory_limit: 1_000_000,
        },
    )
    .unwrap();

    let got = out.lock().unwrap().clone();
    let expected: Vec<_> = (0..20u64).map(|i| (i + 3) * 2 + 1).collect();
    assert_eq!(got, expected);
}

#[test]
fn test_five_stage_pipeline_correctness() {
    // 5 stages of +1, applied to 0..50. Expected: i + 5.
    let out = Arc::new(Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out: out.clone() });
    let source = Box::new(RangeSource { n: 50 });

    let stages: Vec<StageEntry> = (0..5)
        .map(|_| {
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 })))
        })
        .collect();

    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig::default(),
    )
    .unwrap();

    let got = out.lock().unwrap().clone();
    let expected: Vec<_> = (0..50u64).map(|i| i + 5).collect();
    assert_eq!(got, expected);
}

#[test]
fn test_pipeline_cancellation_under_load() {
    let out = Arc::new(Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out: out.clone() });
    let source = Box::new(RangeSource { n: 1_000_000 });
    let counter = Arc::new(AtomicU64::new(0));

    let stages: Vec<StageEntry> = vec![
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(CountingStage {
            counter: counter.clone(),
        }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
    ];

    let cancel = CancelToken::new();
    let cancel_for_timer = cancel.clone();
    std::thread::spawn(move || {
        std::thread::sleep(std::time::Duration::from_millis(50));
        cancel_for_timer.cancel();
    });

    run_multi_stage::<u64, u64>(source, stages, sink, &cancel, PipelineConfig::default()).unwrap();

    // We should have cancelled well before all 1M items were processed.
    let processed = counter.load(Ordering::Relaxed);
    assert!(
        processed < 1_000_000,
        "expected cancellation to stop pipeline early, processed {processed}"
    );
}

#[test]
fn test_pipeline_with_tight_queue_capacity() {
    // Queue capacity of 2 forces lots of held-item activity.
    let out = Arc::new(Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out: out.clone() });
    let source = Box::new(RangeSource { n: 100 });

    let stages: Vec<StageEntry> = vec![
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 10 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(MulStage { n: 3 }))),
    ];

    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads: 1,
            queue_capacity: 2,
            queue_memory_limit: 128,
            global_memory_limit: 10_000,
        },
    )
    .unwrap();

    let got = out.lock().unwrap().clone();
    let expected: Vec<_> = (0..100u64).map(|i| (i + 10) * 3).collect();
    assert_eq!(got, expected);
}

/// Verify the counting-adapter-based completion check accepts a
/// well-behaved pipeline where N items go in and N items come out.
#[test]
fn test_pipeline_completion_counts_match() {
    let out = Arc::new(Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out: out.clone() });
    let source = Box::new(RangeSource { n: 42 });

    let stages: Vec<StageEntry> = vec![
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
    ];

    // Ok => produced == consumed was validated internally.
    run_multi_stage::<u64, u64>(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig::default(),
    )
    .unwrap();

    assert_eq!(out.lock().unwrap().len(), 42);
}

use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn pipeline_output_invariant(
        n in 1u64..200,
        num_stages in 1usize..6,
        threads in 1usize..5,
        queue_cap in 2usize..32,
    ) {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out: out.clone() });
        let source = Box::new(RangeSource { n });

        let stages: Vec<StageEntry> = (0..num_stages)
            .map(|_| StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))))
            .collect();

        run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: threads,
                queue_capacity: queue_cap,
                queue_memory_limit: 4096,
                global_memory_limit: 1_000_000,
            },
        ).unwrap();

        let got = out.lock().unwrap().clone();
        let expected: Vec<_> = (0..n).map(|i| i + num_stages as u64).collect();
        prop_assert_eq!(got, expected);
    }
}

//! Integration tests for the ported `Scheduler` (plan 2c-3a).
//!
//! Exercises the combined pool+scheduler behavior end-to-end; unit tests
//! covering each `Scheduler` method live alongside the scheduler module.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use fgumi_lib::runall::engine::cancel::CancelToken;
use fgumi_lib::runall::engine::driver::{PipelineConfig, StageEntry, run_multi_stage};
use fgumi_lib::runall::engine::stage::{Parallelism, Stage};
use fgumi_lib::runall::engine::stage_source::{ErasedStageSource, StageSource};

use super::mock_stages::{CollectSink, RangeSource};

/// A Parallel stage that burns CPU per item to simulate abundant work.
#[derive(Clone)]
struct BusyParallel;

impl Stage for BusyParallel {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        let mut h: u64 = input;
        for _ in 0..5_000 {
            h = h.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
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
        "BusyParallel"
    }
}

/// A Sequential stage that counts how many items it processed.
#[derive(Clone)]
struct SequentialCounter {
    counter: Arc<AtomicUsize>,
}

impl Stage for SequentialCounter {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        self.counter.fetch_add(1, Ordering::SeqCst);
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
        "SequentialCounter"
    }
}

#[test]
fn test_sequential_stage_not_starved_under_load() {
    // Pipeline: BusyParallel -> SequentialCounter -> BusyParallel.
    // The Sequential stage is in the middle, where it would be starved if
    // no worker owned it. With exclusive ownership (task 5), worker 0
    // owns the Sequential stage and keeps it fed, so all N items arrive.
    const N: u64 = 500;
    const EXPECTED_COUNT: usize = 500;

    let seq_counter = Arc::new(AtomicUsize::new(0));

    let out = Arc::new(std::sync::Mutex::new(Vec::new()));
    let sink = Box::new(CollectSink { out: out.clone() });
    let source = Box::new(RangeSource { n: N });

    let stages: Vec<StageEntry> = vec![
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(BusyParallel))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(SequentialCounter {
            counter: seq_counter.clone(),
        }))),
        StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(BusyParallel))),
    ];

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

    // Every input must have transited the Sequential stage exactly once.
    assert_eq!(
        seq_counter.load(Ordering::SeqCst),
        EXPECTED_COUNT,
        "Sequential stage did not process all items — starvation detected",
    );
    // And every input must have reached the sink.
    let results = out.lock().unwrap();
    assert_eq!(results.len(), EXPECTED_COUNT);
}

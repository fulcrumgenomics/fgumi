//! Error-path integration tests for the `pipeline` engine.
//!
//! Cover panics, source errors, stage errors, and watchdog/cancel interaction
//! to ensure the engine exits cleanly (no hangs) and propagates errors
//! predictably.

use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use fgumi_lib::runall::engine::cancel::CancelToken;
use fgumi_lib::runall::engine::deadlock::WatchdogConfig;
use fgumi_lib::runall::engine::driver::{
    PipelineConfig, StageEntry, run_multi_stage, run_multi_stage_with_watchdog, run_single_stage,
};
use fgumi_lib::runall::engine::sink::InputQueue;
use fgumi_lib::runall::engine::source::{OutputQueue, Source};
use fgumi_lib::runall::engine::special_stage::{SpecialStage, TypedSpecialStage};
use fgumi_lib::runall::engine::stage::{Parallelism, Stage};
use fgumi_lib::runall::engine::stage_source::{ErasedStageSource, StageSource};

use super::mock_stages::{AddStage, CollectSink, RangeSource};

/// Hard upper bound on how long a pipeline test is allowed to take. Any run
/// that exceeds this probably indicates a deadlock regression.
const PIPELINE_TEST_TIMEOUT: Duration = Duration::from_secs(2);

/// Run `f` in a helper thread and fail the test if it doesn't finish in time.
fn run_with_timeout<F>(f: F)
where
    F: FnOnce() + Send + 'static,
{
    let (done_tx, done_rx) = std::sync::mpsc::channel::<()>();
    let handle = std::thread::spawn(move || {
        f();
        let _ = done_tx.send(());
    });
    let start = Instant::now();
    if done_rx.recv_timeout(PIPELINE_TEST_TIMEOUT).is_ok() {
        handle.join().expect("worker thread panicked after signalling done");
    } else {
        panic!(
            "pipeline did not exit within {:?} (elapsed: {:?})",
            PIPELINE_TEST_TIMEOUT,
            start.elapsed(),
        );
    }
}

// ============================================================================
// T1: Worker panic mid-process
// ============================================================================

/// Stage that panics on its first `process` call.
#[derive(Clone)]
pub struct PanickingStage;

impl Stage for PanickingStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, _input: u64, _out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        panic!("intentional test panic in PanickingStage::process");
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, _output: &u64) -> usize {
        std::mem::size_of::<u64>()
    }

    fn name(&self) -> &'static str {
        "PanickingStage"
    }
}

#[test]
fn test_single_stage_worker_panic_returns_err_without_deadlock() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(RangeSource { n: 100 });

        let result = run_single_stage(
            source,
            || PanickingStage,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 1,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 100_000,
            },
        );

        assert!(result.is_err(), "expected Err from panicking worker, got {result:?}");
    });
}

#[test]
fn test_multi_stage_worker_panic_returns_err_without_deadlock() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(RangeSource { n: 100 });

        let stages: Vec<StageEntry> = vec![
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(PanickingStage))),
        ];

        let result = run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 1,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 100_000,
            },
        );

        assert!(result.is_err(), "expected Err from panicking worker, got {result:?}");
    });
}

/// Verify the error message contains the original panic string (extracted
/// from the `Box<dyn Any>` payload) rather than a debug stringification of
/// the box itself. Without the F22 fix, the error reads
/// "thread panicked: Any { .. }" instead of the actual `panic!()` string.
#[test]
fn test_multi_stage_panic_message_is_extracted() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(RangeSource { n: 100 });

        let stages: Vec<StageEntry> = vec![StageEntry::Pool(ErasedStageSource::from_source(
            StageSource::Clone(PanickingStage),
        ))];

        let result = run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 1,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 100_000,
            },
        );

        assert!(result.is_err());
        // The driver adds a `stage N failed` context around the worker error.
        // The anyhow alternate-format `{:#}` walks the error chain and
        // surfaces the underlying panic message extracted by
        // `extract_panic_message`.
        let err = result.unwrap_err();
        let chain = format!("{err:#}");
        assert!(
            chain.contains("intentional test panic in PanickingStage::process"),
            "expected original panic message in error chain, got: {chain}"
        );
    });
}

/// Verify the error also names the offending stage so operators can tell
/// which stage panicked without grepping logs. The pool worker catches the
/// panic at the `process_any` call site and prefixes the error with
/// `stage '<name>' panicked:` before bubbling it through the stage-error
/// path.
#[test]
fn test_multi_stage_panic_error_names_the_stage() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(RangeSource { n: 100 });

        let stages: Vec<StageEntry> = vec![
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(PanickingStage))),
        ];

        let result = run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 2,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 100_000,
            },
        );

        assert!(result.is_err(), "expected Err from panicking stage, got {result:?}");
        let err = result.unwrap_err();
        let chain = format!("{err:#}");
        assert!(
            chain.contains("PanickingStage"),
            "expected stage name 'PanickingStage' in error chain, got: {chain}"
        );
        assert!(chain.contains("panicked"), "expected 'panicked' in error chain, got: {chain}");
        assert!(
            chain.contains("intentional test panic in PanickingStage::process"),
            "expected original panic message in error chain, got: {chain}"
        );
    });
}

/// Exercise a Sequential stage that panics. The pool locks the stage in a
/// mutex and calls `process_any`; the panic path differs from the Parallel
/// branch, so we cover both.
#[test]
fn test_multi_stage_sequential_panic_names_the_stage() {
    #[derive(Clone)]
    struct SequentialPanickingStage;
    impl Stage for SequentialPanickingStage {
        type Input = u64;
        type Output = u64;
        fn process(&mut self, _input: u64, _out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
            panic!("intentional sequential panic");
        }
        fn parallelism(&self) -> Parallelism {
            Parallelism::Sequential
        }
        fn output_memory_estimate(&self, _output: &u64) -> usize {
            std::mem::size_of::<u64>()
        }
        fn name(&self) -> &'static str {
            "SequentialPanickingStage"
        }
    }

    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(RangeSource { n: 100 });

        let stages: Vec<StageEntry> = vec![StageEntry::Pool(ErasedStageSource::from_source(
            StageSource::Clone(SequentialPanickingStage),
        ))];

        let result = run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 2,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 100_000,
            },
        );

        assert!(result.is_err(), "expected Err from panicking sequential stage, got {result:?}");
        let err = result.unwrap_err();
        let chain = format!("{err:#}");
        assert!(
            chain.contains("SequentialPanickingStage"),
            "expected stage name in error chain, got: {chain}"
        );
        assert!(
            chain.contains("intentional sequential panic"),
            "expected original panic message in error chain, got: {chain}"
        );
    });
}

/// Exercise a self-threaded (`SpecialStage`) panic. The driver spawns a
/// dedicated thread for each special stage and wraps the body in
/// `run_and_catch_panic_named`, which prefixes the error with the stage
/// name.
#[test]
fn test_special_stage_panic_names_the_stage() {
    struct PanickingSpecial;
    impl SpecialStage for PanickingSpecial {
        type Input = u64;
        type Output = u64;
        fn run(
            self: Box<Self>,
            input: Box<dyn InputQueue<u64>>,
            output: Box<dyn OutputQueue<u64>>,
            _cancel: CancelToken,
        ) -> anyhow::Result<()> {
            // Consume one item so the upstream source can make forward
            // progress, then intentionally panic. Close defensively on the
            // error path so the driver's joins don't block.
            let _ = input.pop();
            output.close();
            panic!("intentional panic in PanickingSpecial::run");
        }
        fn name(&self) -> &'static str {
            "PanickingSpecial"
        }
    }

    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(RangeSource { n: 100 });

        let stages: Vec<StageEntry> =
            vec![StageEntry::Special(Box::new(TypedSpecialStage::new(PanickingSpecial)))];

        let result = run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 1,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 100_000,
            },
        );

        assert!(result.is_err(), "expected Err from panicking special stage, got {result:?}");
        let err = result.unwrap_err();
        let chain = format!("{err:#}");
        assert!(
            chain.contains("PanickingSpecial"),
            "expected special-stage name 'PanickingSpecial' in error chain, got: {chain}"
        );
        assert!(
            chain.contains("intentional panic in PanickingSpecial::run"),
            "expected original panic message in error chain, got: {chain}"
        );
    });
}

// ============================================================================
// T2: Source returns Err
// ============================================================================

/// Source whose `run` immediately returns an error.
pub struct FailingSource;

impl Source for FailingSource {
    type Output = u64;

    fn run(
        self: Box<Self>,
        output: Box<dyn OutputQueue<Self::Output>>,
        _cancel: CancelToken,
    ) -> anyhow::Result<()> {
        output.close();
        anyhow::bail!("boom");
    }

    fn name(&self) -> &'static str {
        "FailingSource"
    }
}

#[test]
fn test_single_stage_source_error_is_propagated() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(FailingSource);

        let result = run_single_stage(
            source,
            || AddStage { n: 1 },
            sink,
            &CancelToken::new(),
            PipelineConfig {
                worker_threads: 1,
                queue_capacity: 4,
                queue_memory_limit: 1024,
                global_memory_limit: 100_000,
            },
        );

        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("boom"), "expected 'boom' in error, got: {msg}");
    });
}

#[test]
fn test_multi_stage_source_error_is_propagated() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(FailingSource);

        let stages: Vec<StageEntry> =
            vec![StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage {
                n: 1,
            })))];

        let result = run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig::default(),
        );

        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("boom"), "expected 'boom' in error, got: {msg}");
    });
}

// ============================================================================
// T3: Stage returns Err in multi-stage path
// ============================================================================

/// Stage that returns an error after N successful items.
#[derive(Clone)]
pub struct FailingStage {
    pub fail_after: u64,
    pub counter: Arc<AtomicU64>,
}

impl Stage for FailingStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        let n = self.counter.fetch_add(1, Ordering::SeqCst);
        if n >= self.fail_after {
            anyhow::bail!("boom");
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
        "FailingStage"
    }
}

#[test]
fn test_multi_stage_stage_error_is_propagated() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out });
        let source = Box::new(RangeSource { n: 100 });
        let counter = Arc::new(AtomicU64::new(0));

        let stages: Vec<StageEntry> = vec![
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage { n: 1 }))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(FailingStage {
                fail_after: 5,
                counter: counter.clone(),
            }))),
        ];

        let result = run_multi_stage::<u64, u64>(
            source,
            stages,
            sink,
            &CancelToken::new(),
            PipelineConfig::default(),
        );

        assert!(result.is_err(), "expected Err from failing stage, got {result:?}");
        let err = result.unwrap_err();
        // The stage error is wrapped with "stage N failed" context; the
        // original "boom" is available via the cause chain.
        let chain: Vec<String> = err.chain().map(ToString::to_string).collect();
        assert!(
            chain.iter().any(|m| m.contains("boom")),
            "expected 'boom' somewhere in error chain, got: {chain:?}",
        );
    });
}

// ============================================================================
// T4: Watchdog with a SUCCESSFUL pipeline
// ============================================================================

#[test]
fn test_watchdog_does_not_cancel_caller_token_on_success() {
    run_with_timeout(|| {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(CollectSink { out: out.clone() });
        let source = Box::new(RangeSource { n: 50 });

        let stages: Vec<StageEntry> =
            vec![StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(AddStage {
                n: 1,
            })))];

        let cancel = CancelToken::new();
        let watchdog_cfg = WatchdogConfig {
            timeout: Duration::from_secs(5),
            poll_interval: Duration::from_millis(20),
        };

        let result = run_multi_stage_with_watchdog::<u64, u64>(
            source,
            stages,
            sink,
            &cancel,
            PipelineConfig::default(),
            Some(watchdog_cfg),
        );

        assert!(result.is_ok(), "expected success, got {result:?}");
        assert!(
            !cancel.is_cancelled(),
            "watchdog wrapper must NOT poison caller's cancel token on success",
        );

        let collected = out.lock().unwrap().clone();
        let expected: Vec<_> = (0..50u64).map(|i| i + 1).collect();
        assert_eq!(collected, expected);
    });
}

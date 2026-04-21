//! # Pipeline Builder
//!
//! Fluent, type-state builder for composing and running pipelines.
//! Wraps the lower-level `run_multi_stage` entry point (see
//! [`crate::runall::engine::driver`]) with a call order that the type
//! system can enforce.
//!
//! ```ignore
//! Pipeline::builder()
//!     .source(my_source)
//!     .stage(stage_a)
//!     .stage(stage_b)
//!     .sink(my_sink)
//!     .threads(8)
//!     .build()
//!     .run(cancel)?;
//! ```
//!
//! ## Type-state call order
//!
//! Each transition consumes the builder and returns a new state. The
//! compiler enforces:
//!
//! 1. `PipelineBuilder::new()` → call `.source(...)` first.
//! 2. [`PipelineBuilderNoStage`] — source set, no stages yet. Must
//!    call `.stage(...)`, `.stage_factory(...)`, or
//!    `.special_stage(...)` at least once. `.sink(...)` is **not**
//!    yet reachable (a pipeline with zero stages is a compile error,
//!    not a runtime error).
//! 3. [`PipelineBuilderWithStage`] — source and ≥1 stage set.
//!    Additional stages may be appended; each stage's `Input` type
//!    must match the prior stage's `Output` (enforced by the type
//!    state). Finalize by calling `.sink(...)` to get a
//!    [`Pipeline`].
//!
//! At any state, `.config(...)` replaces the full configuration and
//! `.threads(N)` overrides `PipelineConfig::worker_threads`.
//!
//! ## Stage variants
//!
//! - `.stage(st)` — requires `Stage: Clone + Send + Sync + 'static`.
//!   The prototype is cloned once per pool worker; worker N holds
//!   `st.clone()`.
//! - `.stage_factory(f)` — for stages whose state cannot be `Clone`
//!   (persistent compressor, `Box<dyn Trait>` fields, etc.). The
//!   factory produces a fresh instance per worker. Type-safe:
//!   `St::Input` must match the current output type.
//! - `.special_stage(s)` — appends a
//!   [`SpecialStage`] that will own a dedicated OS thread at run time.
//!
//! ## Runtime topology
//!
//! The builder and any direct `run_multi_stage` caller share the same
//! runtime: source → pool-driven stages + dedicated-thread special
//! stages → sink, with `StageQueue` transitions between stages,
//! `ReorderBuffer`-based order restoration where needed, and shared
//! `CancelToken` / `MemoryTracker`. See the [driver](super::driver)
//! module docs for the full topology, cancellation, and memory-gating
//! overview.

use super::cancel::CancelToken;
use super::driver::{ErasedStage, PipelineConfig, StageEntry, TypedStage, run_multi_stage};
use super::sink::Sink;
use super::source::Source;
use super::special_stage::{SpecialStage, TypedSpecialStage};
use super::stage::Stage;
use super::stage_source::{ErasedStageSource, StageSource};

/// A configured, ready-to-run pipeline.
pub struct Pipeline<SourceOut: Send + 'static, SinkIn: Send + 'static> {
    source: Box<dyn Source<Output = SourceOut>>,
    stages: Vec<StageEntry>,
    sink: Box<dyn Sink<Input = SinkIn>>,
    config: PipelineConfig,
}

impl<SourceOut: Send + 'static, SinkIn: Send + 'static> Pipeline<SourceOut, SinkIn> {
    /// Run the configured pipeline to completion.
    ///
    /// # Errors
    ///
    /// Returns an error if any pipeline thread fails or panics, or if there are
    /// no stages configured.
    pub fn run(self, cancel: &CancelToken) -> anyhow::Result<()> {
        run_multi_stage(self.source, self.stages, self.sink, cancel, self.config)
    }
}

/// Builder state before the source is set.
pub struct PipelineBuilder;

/// Builder state with source set but NO stages yet.
///
/// Exposes `.stage(...)`, `.config(...)`, and `.threads(...)` — but NOT
/// `.sink(...)`. This makes a pipeline-with-zero-stages a compile error
/// rather than a runtime error (`run_multi_stage` requires at least one
/// stage).
pub struct PipelineBuilderNoStage<SourceOut: Send + 'static> {
    source: Box<dyn Source<Output = SourceOut>>,
    config: PipelineConfig,
}

/// Builder state with source AND at least one stage set. Adds `.sink(...)`
/// to complete the pipeline.
pub struct PipelineBuilderWithStage<SourceOut: Send + 'static, CurrentOut: Send + 'static> {
    source: Box<dyn Source<Output = SourceOut>>,
    stages: Vec<StageEntry>,
    config: PipelineConfig,
    _marker: std::marker::PhantomData<CurrentOut>,
}

impl PipelineBuilder {
    /// Create a new empty builder.
    #[must_use]
    pub fn new() -> Self {
        Self
    }

    /// Set the pipeline's source.
    ///
    /// Returns a builder that only exposes `.stage(...)` — pipelines must
    /// contain at least one stage, enforced by the type state.
    #[must_use]
    pub fn source<SourceOut: Send + 'static>(
        self,
        source: Box<dyn Source<Output = SourceOut>>,
    ) -> PipelineBuilderNoStage<SourceOut> {
        PipelineBuilderNoStage { source, config: PipelineConfig::default() }
    }
}

impl Default for PipelineBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl<SourceOut> PipelineBuilderNoStage<SourceOut>
where
    SourceOut: Send + 'static,
{
    /// Add the first stage. The stage's `Input` type must match the source's
    /// `Output` type. Transitions the builder into a state where subsequent
    /// `.stage(...)` calls are allowed and `.sink(...)` becomes reachable.
    ///
    /// The stage prototype is cloned once per pool worker. For stages whose
    /// state cannot be `Clone`, use [`Self::stage_factory`] instead.
    #[must_use]
    pub fn stage<St>(self, stage: St) -> PipelineBuilderWithStage<SourceOut, St::Output>
    where
        St: Stage<Input = SourceOut> + Clone + Send + Sync + 'static,
    {
        PipelineBuilderWithStage {
            source: self.source,
            stages: vec![StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(
                stage,
            )))],
            config: self.config,
            _marker: std::marker::PhantomData,
        }
    }

    /// Add the first stage using a factory that produces fresh instances per
    /// worker. Use this for stages that own non-`Clone` state (for example,
    /// a stage holding a persistent compressor or a `Box<dyn Trait>` field).
    ///
    /// Type-safe: the compiler enforces `St::Input == SourceOut` at the call
    /// site. No `Clone` bound — the factory itself produces fresh instances.
    #[must_use]
    pub fn stage_factory<St, F>(self, factory: F) -> PipelineBuilderWithStage<SourceOut, St::Output>
    where
        St: Stage<Input = SourceOut> + Send + 'static,
        F: Fn() -> St + Send + Sync + 'static,
    {
        let probe = factory();
        let parallelism = probe.parallelism();
        let name = probe.name();
        drop(probe);
        PipelineBuilderWithStage {
            source: self.source,
            stages: vec![StageEntry::Pool(ErasedStageSource::from_erased_factory(
                move || Box::new(TypedStage::new(factory())) as Box<dyn ErasedStage>,
                parallelism,
                name,
            ))],
            config: self.config,
            _marker: std::marker::PhantomData,
        }
    }

    /// Add the first stage as a [`SpecialStage`].
    ///
    /// The stage's `Input` type must match the source's `Output` type.
    /// Transitions the builder to `PipelineBuilderWithStage`.
    #[must_use]
    pub fn special_stage<S>(self, stage: S) -> PipelineBuilderWithStage<SourceOut, S::Output>
    where
        S: SpecialStage<Input = SourceOut> + 'static,
    {
        PipelineBuilderWithStage {
            source: self.source,
            stages: vec![StageEntry::Special(Box::new(TypedSpecialStage::new(stage)))],
            config: self.config,
            _marker: std::marker::PhantomData,
        }
    }

    /// Replace the entire pipeline configuration.
    #[must_use]
    pub fn config(mut self, config: PipelineConfig) -> Self {
        self.config = config;
        self
    }

    /// Override the worker thread count.
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.config.worker_threads = threads;
        self
    }
}

impl<SourceOut, CurrentOut> PipelineBuilderWithStage<SourceOut, CurrentOut>
where
    SourceOut: Send + 'static,
    CurrentOut: Send + 'static,
{
    /// Append another stage. The stage's `Input` type must match the current
    /// output type.
    ///
    /// The stage prototype is cloned once per pool worker. For stages whose
    /// state cannot be `Clone`, use [`Self::stage_factory`] instead.
    #[must_use]
    pub fn stage<St>(mut self, stage: St) -> PipelineBuilderWithStage<SourceOut, St::Output>
    where
        St: Stage<Input = CurrentOut> + Clone + Send + Sync + 'static,
    {
        self.stages
            .push(StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(stage))));
        PipelineBuilderWithStage {
            source: self.source,
            stages: self.stages,
            config: self.config,
            _marker: std::marker::PhantomData,
        }
    }

    /// Append another stage using a factory that produces fresh instances per
    /// worker. Use this for stages that own non-`Clone` state.
    ///
    /// Type-safe: the compiler enforces `St::Input == CurrentOut` at the call
    /// site. No `Clone` bound — the factory itself produces fresh instances.
    #[must_use]
    pub fn stage_factory<St, F>(
        mut self,
        factory: F,
    ) -> PipelineBuilderWithStage<SourceOut, St::Output>
    where
        St: Stage<Input = CurrentOut> + Send + 'static,
        F: Fn() -> St + Send + Sync + 'static,
    {
        let probe = factory();
        let parallelism = probe.parallelism();
        let name = probe.name();
        drop(probe);
        self.stages.push(StageEntry::Pool(ErasedStageSource::from_erased_factory(
            move || Box::new(TypedStage::new(factory())) as Box<dyn ErasedStage>,
            parallelism,
            name,
        )));
        PipelineBuilderWithStage {
            source: self.source,
            stages: self.stages,
            config: self.config,
            _marker: std::marker::PhantomData,
        }
    }

    /// Append a [`SpecialStage`]. The stage's `Input` type must match the
    /// current output type.
    #[must_use]
    pub fn special_stage<S>(mut self, stage: S) -> PipelineBuilderWithStage<SourceOut, S::Output>
    where
        S: SpecialStage<Input = CurrentOut> + 'static,
    {
        self.stages.push(StageEntry::Special(Box::new(TypedSpecialStage::new(stage))));
        PipelineBuilderWithStage {
            source: self.source,
            stages: self.stages,
            config: self.config,
            _marker: std::marker::PhantomData,
        }
    }

    /// Replace the entire pipeline configuration.
    #[must_use]
    pub fn config(mut self, config: PipelineConfig) -> Self {
        self.config = config;
        self
    }

    /// Override the worker thread count.
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.config.worker_threads = threads;
        self
    }

    /// Finalize the pipeline by setting the sink.
    #[must_use]
    pub fn sink(self, sink: Box<dyn Sink<Input = CurrentOut>>) -> Pipeline<SourceOut, CurrentOut> {
        Pipeline { source: self.source, stages: self.stages, sink, config: self.config }
    }
}

impl Pipeline<(), ()> {
    /// Begin building a new pipeline.
    #[must_use]
    pub fn builder() -> PipelineBuilder {
        PipelineBuilder::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::cancel::CancelToken;
    use crate::runall::engine::sink::InputQueue;
    use crate::runall::engine::source::OutputQueue;
    use crate::runall::engine::stage::{Parallelism, SequencedItem};
    use std::sync::{Arc, Mutex};

    struct Numbers {
        n: u64,
    }
    impl Source for Numbers {
        type Output = u64;
        fn run(
            self: Box<Self>,
            output: Box<dyn OutputQueue<u64>>,
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
            "Numbers"
        }
    }

    #[derive(Clone)]
    struct Plus1;
    impl Stage for Plus1 {
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
            "Plus1"
        }
    }

    struct Collect {
        out: Arc<Mutex<Vec<u64>>>,
    }
    impl Sink for Collect {
        type Input = u64;
        fn run(
            self: Box<Self>,
            input: Box<dyn InputQueue<u64>>,
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
            buffer.sort_by_key(|(s, _)| *s);
            let mut o = self.out.lock().unwrap();
            for (_, v) in buffer {
                o.push(v);
            }
            Ok(())
        }
        fn name(&self) -> &'static str {
            "Collect"
        }
    }

    #[test]
    fn test_builder_composes_pipeline() {
        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(Collect { out: out.clone() });

        let pipeline = PipelineBuilder::new()
            .source(Box::new(Numbers { n: 10 }))
            .stage(Plus1)
            .stage(Plus1)
            .threads(2)
            .sink(sink);

        pipeline.run(&CancelToken::new()).unwrap();

        let expected: Vec<_> = (0..10u64).map(|i| i + 2).collect();
        assert_eq!(*out.lock().unwrap(), expected);
    }

    /// [`SpecialStage`] integration: Source -> pool stage -> special stage (passthrough) -> pool stage -> Sink.
    #[test]
    fn test_special_stage_engine_integration() {
        use crate::runall::engine::special_stage::SpecialStage;

        /// Passthrough special stage: copies input to output unchanged.
        struct Passthrough;

        impl SpecialStage for Passthrough {
            type Input = u64;
            type Output = u64;

            fn run(
                self: Box<Self>,
                input: Box<dyn InputQueue<u64>>,
                output: Box<dyn OutputQueue<u64>>,
                cancel: CancelToken,
            ) -> anyhow::Result<()> {
                loop {
                    if cancel.is_cancelled() {
                        break;
                    }
                    if let Some(item) = input.pop() {
                        let out_item =
                            SequencedItem::new(item.seq, item.item, item.memory_estimate);
                        if output.push_until_cancelled(out_item, &cancel).is_err() {
                            break;
                        }
                    } else if input.is_drained() {
                        break;
                    } else {
                        std::thread::yield_now();
                    }
                }
                output.close();
                Ok(())
            }

            fn name(&self) -> &'static str {
                "Passthrough"
            }
        }

        let out = Arc::new(Mutex::new(Vec::new()));
        let sink = Box::new(Collect { out: out.clone() });

        let pipeline = PipelineBuilder::new()
            .source(Box::new(Numbers { n: 20 }))
            .stage(Plus1) // pool stage: +1
            .special_stage(Passthrough) // special stage: passthrough
            .stage(Plus1) // pool stage: +1
            .threads(2)
            .sink(sink);

        pipeline.run(&CancelToken::new()).unwrap();

        let mut expected: Vec<_> = (0..20u64).map(|i| i + 2).collect();
        let mut result = out.lock().unwrap().clone();
        // Pool stages may reorder; sort both for comparison.
        expected.sort_unstable();
        result.sort_unstable();
        assert_eq!(result, expected);
    }
}

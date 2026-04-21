//! Policy for how the work-stealing pool constructs per-worker stage instances.
//!
//! Parallel stages are executed concurrently by multiple workers. Each worker
//! needs its own stage instance because `Stage::process` takes `&mut self`.
//!
//! `StageSource<S>` captures the two supported construction policies:
//!
//! - `Clone(S)`: stage is `Clone`. Each worker clones the original once at
//!   pool start. Simplest; use for stateless or Clone-safe stages.
//! - `Factory(Arc<Fn() -> S>)`: stage has non-clone state (for example a
//!   `Box<dyn ConsensusCaller>` that cannot be cloned). The factory produces
//!   fresh instances per worker.
//!
//! For Sequential stages, the pool uses ONE instance wrapped in a `Mutex`
//! regardless of policy (since only one worker runs it at a time).
//! For Barrier and Special stages, stages are handled separately (not in the
//! pool).

use std::sync::Arc;

use super::driver::{ErasedStage, TypedStage};
use super::stage::{Parallelism, Stage};

/// How the pool obtains per-worker stage instances for a given logical stage.
pub enum StageSource<S: Stage + 'static> {
    /// A `Clone`-based stage. Cloned once per worker at pool start.
    Clone(S),
    /// A factory that produces fresh stage instances on demand.
    Factory(Arc<dyn Fn() -> S + Send + Sync>),
}

impl<S: Stage + Clone + 'static> From<S> for StageSource<S> {
    fn from(stage: S) -> Self {
        StageSource::Clone(stage)
    }
}

impl<S: Stage + 'static> StageSource<S> {
    /// Build a factory-based source from a closure.
    pub fn factory(f: impl Fn() -> S + Send + Sync + 'static) -> Self {
        StageSource::Factory(Arc::new(f))
    }

    /// Build a per-worker instance.
    ///
    /// Requires `S: Clone` for the `Clone` variant. Factories do not need
    /// `Clone` because each call produces a fresh instance.
    pub fn build(&self) -> S
    where
        S: Clone,
    {
        match self {
            StageSource::Clone(s) => s.clone(),
            StageSource::Factory(f) => f(),
        }
    }

    /// Build a per-worker instance when `S` may not be `Clone`.
    ///
    /// Returns `None` for the `Clone` variant (caller error — the `Clone`
    /// variant requires `S: Clone`, which this signature cannot enforce).
    pub fn build_unchecked(&self) -> Option<S> {
        match self {
            StageSource::Factory(f) => Some(f()),
            StageSource::Clone(_) => None,
        }
    }
}

/// Type-erased wrapper around a [`StageSource`]. The engine stores a
/// `Vec<ErasedStageSource>` and invokes [`ErasedStageSource::build_erased`] to
/// obtain per-worker `Box<dyn ErasedStage>` instances.
pub struct ErasedStageSource {
    /// Construction function. Called per worker for parallel stages; called
    /// once for Sequential.
    builder: Arc<dyn Fn() -> Box<dyn ErasedStage> + Send + Sync>,
    /// Parallelism declared by the underlying stage.
    parallelism: Parallelism,
    /// Human-readable name for logs and diagnostics.
    name: &'static str,
}

impl ErasedStageSource {
    /// Construct an erased source from a typed [`StageSource<S>`]. The
    /// underlying stage must be `Clone + Send + Sync` so that both `Clone`
    /// and `Factory` variants can be shared safely across pool workers and
    /// uniformly hand out fresh instances.
    #[must_use]
    pub fn from_source<S>(source: StageSource<S>) -> Self
    where
        S: Stage + Clone + Send + Sync + 'static,
    {
        match source {
            StageSource::Clone(prototype) => {
                let parallelism = prototype.parallelism();
                let name = prototype.name();
                // Capture the prototype by move into the closure; clone on
                // each call. Avoids any `Arc<StageSource<S>>` Sync bounds.
                let builder: Arc<dyn Fn() -> Box<dyn ErasedStage> + Send + Sync> =
                    Arc::new(move || {
                        let instance = prototype.clone();
                        Box::new(TypedStage::new(instance)) as Box<dyn ErasedStage>
                    });
                Self { builder, parallelism, name }
            }
            StageSource::Factory(factory) => {
                // Factories produce fresh Parallel stages. The factory is
                // already Send + Sync so we just wrap invocation.
                let builder: Arc<dyn Fn() -> Box<dyn ErasedStage> + Send + Sync> =
                    Arc::new(move || {
                        let instance = factory();
                        Box::new(TypedStage::new(instance)) as Box<dyn ErasedStage>
                    });
                Self { builder, parallelism: Parallelism::Parallel, name: "Factory" }
            }
        }
    }

    /// Construct from an arbitrary erased-stage factory. Use this for stages
    /// whose state cannot be expressed as `S: Stage + Clone` (for example,
    /// stages that own a `Box<dyn Trait>` field).
    #[must_use]
    pub fn from_erased_factory(
        factory: impl Fn() -> Box<dyn ErasedStage> + Send + Sync + 'static,
        parallelism: Parallelism,
        name: &'static str,
    ) -> Self {
        Self { builder: Arc::new(factory), parallelism, name }
    }

    /// Produce a fresh erased stage instance.
    #[must_use]
    pub fn build_erased(&self) -> Box<dyn ErasedStage> {
        (self.builder)()
    }

    /// Declared parallelism of this stage.
    #[must_use]
    pub fn parallelism(&self) -> Parallelism {
        self.parallelism
    }

    /// Human-readable stage name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        self.name
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone)]
    struct CloneableStage {
        n: u32,
    }

    impl Stage for CloneableStage {
        type Input = u64;
        type Output = u64;

        fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
            self.n += 1;
            out(input + u64::from(self.n));
            Ok(())
        }

        fn parallelism(&self) -> Parallelism {
            Parallelism::Parallel
        }

        fn output_memory_estimate(&self, _o: &u64) -> usize {
            8
        }

        fn name(&self) -> &'static str {
            "CloneableStage"
        }
    }

    #[test]
    fn test_clone_source_builds_via_clone() {
        let src = StageSource::Clone(CloneableStage { n: 5 });
        let s1 = src.build();
        assert_eq!(s1.n, 5);
    }

    #[test]
    fn test_factory_source_builds_via_factory() {
        let src = StageSource::factory(|| CloneableStage { n: 42 });
        let s1 = src.build();
        assert_eq!(s1.n, 42);
    }

    #[test]
    fn test_from_impl_wraps_clone_stages() {
        let s = CloneableStage { n: 1 };
        let src: StageSource<CloneableStage> = s.into();
        assert!(matches!(src, StageSource::Clone(_)));
    }

    #[test]
    fn test_erased_source_from_clone_stage() {
        let src = ErasedStageSource::from_source(StageSource::Clone(CloneableStage { n: 42 }));
        assert_eq!(src.parallelism(), Parallelism::Parallel);
        let _stage1 = src.build_erased();
        let _stage2 = src.build_erased(); // multiple builds OK
    }

    #[test]
    fn test_erased_source_from_factory() {
        let counter = std::sync::Arc::new(std::sync::atomic::AtomicU32::new(0));
        let c = counter.clone();
        let src =
            ErasedStageSource::from_source::<CloneableStage>(StageSource::factory(move || {
                let id = c.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                CloneableStage { n: id }
            }));
        let _s1 = src.build_erased();
        let _s2 = src.build_erased();
        assert_eq!(counter.load(std::sync::atomic::Ordering::SeqCst), 2);
    }

    /// Sequential stages carry their parallelism declaration through erasure
    /// so the pool can decide to install a Mutex instead of per-worker clones.
    #[test]
    fn test_erased_source_preserves_sequential_parallelism() {
        #[derive(Clone)]
        struct SeqStage;
        impl Stage for SeqStage {
            type Input = u64;
            type Output = u64;
            fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
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
                "SeqStage"
            }
        }
        let src = ErasedStageSource::from_source(StageSource::Clone(SeqStage));
        assert_eq!(src.parallelism(), Parallelism::Sequential);
        assert_eq!(src.name(), "SeqStage");
    }

    /// The raw erased-factory constructor preserves its caller-provided
    /// parallelism and name so the engine can build heterogeneous stage
    /// chains without knowing each stage's concrete Rust type.
    #[test]
    fn test_from_erased_factory_preserves_declared_parallelism_and_name() {
        let src = ErasedStageSource::from_erased_factory(
            || Box::new(TypedStage::new(CloneableStage { n: 7 })) as Box<dyn ErasedStage>,
            Parallelism::Parallel,
            "erased_test",
        );
        assert_eq!(src.parallelism(), Parallelism::Parallel);
        assert_eq!(src.name(), "erased_test");
        let _s = src.build_erased();
    }
}

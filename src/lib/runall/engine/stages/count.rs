//! Dummy pass-through counter stage. Stand-in for real `ExtractStage`
//! during plan 2c-1. Increments a shared counter by the number of
//! templates in each batch.
//!
//! ## Input
//!
//! [`TemplateBatch`] — a batch of FASTQ templates.
//!
//! ## Output
//!
//! `()` — no downstream payload; this stage is purely a side-effecting
//! sink-adjacent probe.
//!
//! ## Ordering guarantees
//!
//! None: counts are aggregated via a relaxed [`AtomicU64`]. The order
//! in which batches increment the counter is not observable in the
//! final value.
//!
//! ## Memory model
//!
//! Zero additional heap allocation per batch.
//!
//! ## Determinism
//!
//! Final counter value is deterministic (addition is associative and
//! commutative). No byte-level output to compare.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;

use crate::runall::engine::fastq_types::TemplateBatch;
use crate::runall::engine::stage::{Parallelism, Stage};

#[derive(Clone)]
pub struct CountStage {
    counter: Arc<AtomicU64>,
}

impl CountStage {
    #[must_use]
    pub fn new(counter: Arc<AtomicU64>) -> Self {
        Self { counter }
    }
}

impl Stage for CountStage {
    type Input = TemplateBatch;
    type Output = ();

    #[tracing::instrument(name = "count", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let n = u64::try_from(input.templates.len()).unwrap_or(u64::MAX);
        self.counter.fetch_add(n, Ordering::Relaxed);
        out(());
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, (): &Self::Output) -> usize {
        0
    }

    fn name(&self) -> &'static str {
        "CountStage"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::fastq_types::{FastqTemplateV2, OwnedFastqRecord, TemplateBatch};

    fn template(name: &[u8]) -> FastqTemplateV2 {
        FastqTemplateV2 {
            records: vec![OwnedFastqRecord {
                name: name.to_vec(),
                sequence: b"A".to_vec(),
                quality: b"I".to_vec(),
            }],
            name: name.to_vec(),
        }
    }

    fn batch(templates: Vec<FastqTemplateV2>) -> TemplateBatch {
        TemplateBatch { templates, ordinal: 0 }
    }

    #[test]
    fn test_count_stage_increments_shared_counter() {
        let counter = Arc::new(AtomicU64::new(0));
        let mut stage = CountStage::new(counter.clone());
        stage
            .process(batch(vec![template(b"r1"), template(b"r2"), template(b"r3")]), &mut |()| {})
            .unwrap();
        assert_eq!(counter.load(Ordering::Relaxed), 3);
        stage.process(batch(vec![template(b"r4"), template(b"r5")]), &mut |()| {}).unwrap();
        assert_eq!(counter.load(Ordering::Relaxed), 5);
    }

    #[test]
    fn test_count_stage_empty_batch() {
        let counter = Arc::new(AtomicU64::new(0));
        let mut stage = CountStage::new(counter.clone());
        stage.process(batch(vec![]), &mut |()| {}).unwrap();
        assert_eq!(counter.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_count_stage_clone_shares_counter() {
        let counter = Arc::new(AtomicU64::new(0));
        let mut stage1 = CountStage::new(counter.clone());
        let mut stage2 = stage1.clone();
        stage1.process(batch(vec![template(b"r1")]), &mut |()| {}).unwrap();
        stage2.process(batch(vec![template(b"r2"), template(b"r3")]), &mut |()| {}).unwrap();
        assert_eq!(counter.load(Ordering::Relaxed), 3);
    }
}

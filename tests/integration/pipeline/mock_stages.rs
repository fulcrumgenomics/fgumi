//! Mock stage implementations for `pipeline` engine tests.
//!
//! These stages operate on `u64` and perform simple arithmetic so we can
//! validate pipeline mechanics without coupling to domain logic.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use fgumi_lib::runall::engine::cancel::CancelToken;
use fgumi_lib::runall::engine::sink::{InputQueue, Sink};
use fgumi_lib::runall::engine::source::{OutputQueue, Source};
use fgumi_lib::runall::engine::stage::{Parallelism, SequencedItem, Stage};

/// Source that emits `n` consecutive u64 values starting from 0.
pub struct RangeSource {
    pub n: u64,
}

impl Source for RangeSource {
    type Output = u64;

    fn run(
        self: Box<Self>,
        output: Box<dyn OutputQueue<Self::Output>>,
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
        "RangeSource"
    }
}

/// Stage that adds a constant to the input.
#[derive(Clone)]
pub struct AddStage {
    pub n: u64,
}

impl Stage for AddStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        out(input + self.n);
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, _output: &u64) -> usize {
        std::mem::size_of::<u64>()
    }

    fn name(&self) -> &'static str {
        "AddStage"
    }
}

/// Stage that multiplies the input by a constant.
#[derive(Clone)]
pub struct MulStage {
    pub n: u64,
}

impl Stage for MulStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        out(input * self.n);
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, _output: &u64) -> usize {
        std::mem::size_of::<u64>()
    }

    fn name(&self) -> &'static str {
        "MulStage"
    }
}

/// Stage that counts items it's seen (for liveness checks).
#[derive(Clone)]
pub struct CountingStage {
    pub counter: Arc<AtomicU64>,
}

impl Stage for CountingStage {
    type Input = u64;
    type Output = u64;

    fn process(&mut self, input: u64, out: &mut dyn FnMut(u64)) -> anyhow::Result<()> {
        self.counter.fetch_add(1, Ordering::Relaxed);
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
        "CountingStage"
    }
}

/// Sink that collects items in sequence-number order.
pub struct CollectSink {
    pub out: Arc<std::sync::Mutex<Vec<u64>>>,
}

impl Sink for CollectSink {
    type Input = u64;

    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
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
        buffer.sort_by_key(|(seq, _)| *seq);
        let mut out = self.out.lock().unwrap();
        for (_, v) in buffer {
            out.push(v);
        }
        Ok(())
    }

    fn name(&self) -> &'static str {
        "CollectSink"
    }
}

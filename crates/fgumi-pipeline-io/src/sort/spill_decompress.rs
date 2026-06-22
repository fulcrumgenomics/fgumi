//! `SortSpillDecompress` — Parallel typed step that reads spill chunk
//! files, decompresses them inline, and pushes the decompressed bytes
//! into per-slot bounded queues on `SortMergeSlot`.

use std::io;
use std::sync::Arc;

use fgumi_sort::{PHASE2_DECOMP_CAP, SortMergeSlot, SpillBlockDecompressor};
use parking_lot::Mutex;

use crate::sort::protocol::{SortPhase1Event, SortPhase2Event};
use fgumi_pipeline_core::{
    Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Max raw blocks read+decompressed by a single `try_fill_some_slot` call.
const MAX_BATCH_PER_CALL: usize = 4;

struct RegisteredSpill {
    slot: Arc<SortMergeSlot>,
}

/// Parallel step that reads + decompresses spill chunk files and
/// pushes results into per-slot queues. Forwards `SortPhase1Event`s
/// verbatim to `SortMerge`.
pub struct SortSpillDecompress {
    registry: Arc<Mutex<Vec<RegisteredSpill>>>,
    block_dec: SpillBlockDecompressor,
    held: HeldSlot<Unpushed<SortPhase2Event>>,
    output_byte_limit: u64,
}

impl SortSpillDecompress {
    /// Construct a fresh step with an empty registry.
    ///
    /// `output_byte_limit` byte-bounds the forwarded-event output queue.
    /// The forwarded `SortPhase2Event::MemoryChunk` variant retains sorted
    /// record chunks, so this queue must budget on bytes (`HeapSize`), not
    /// event count, to keep retained memory a function of configuration.
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            registry: Arc::new(Mutex::new(Vec::new())),
            block_dec: SpillBlockDecompressor::new(),
            held: HeldSlot::new(),
            output_byte_limit,
        }
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        let Some(unpushed) = self.held.take() else {
            return true;
        };
        match ctx.outputs.retry(unpushed) {
            Ok(()) => true,
            Err(again) => {
                self.held.put(again);
                false
            }
        }
    }

    fn push_or_hold(&mut self, ctx: &mut StepCtx<'_, Self>, event: SortPhase2Event) -> bool {
        match ctx.outputs.push(event) {
            Ok(()) => true,
            Err(unpushed) => {
                self.held.put(unpushed);
                false
            }
        }
    }

    fn snapshot_registry(&self) -> Vec<Arc<SortMergeSlot>> {
        let registry = self.registry.lock();
        registry.iter().map(|e| Arc::clone(&e.slot)).collect()
    }

    fn try_fill_some_slot(&mut self) -> io::Result<bool> {
        for slot in self.snapshot_registry() {
            if slot.queue_eof.load(std::sync::atomic::Ordering::Acquire) {
                continue;
            }

            let Ok(mut reader_guard) = slot.reader.try_lock() else {
                continue;
            };

            let m = {
                let dec = slot.decompressed.lock().expect("decompressed mutex poisoned");
                PHASE2_DECOMP_CAP.saturating_sub(dec.len())
            };
            if m == 0 {
                continue;
            }
            let m = m.min(MAX_BATCH_PER_CALL);

            let decompressed_batch =
                match self.block_dec.read_blocks(&mut reader_guard.inner, slot.codec, m) {
                    Ok(b) => b,
                    Err(e) => {
                        {
                            let _g = slot.decompressed.lock().expect("decompressed mutex poisoned");
                            slot.decomp_error.store(true, std::sync::atomic::Ordering::Release);
                            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
                        }
                        drop(reader_guard);
                        return Err(e);
                    }
                };
            let got = decompressed_batch.len();
            let hit_eof = got < m;

            if got == 0 {
                {
                    let _g = slot.decompressed.lock().expect("decompressed mutex poisoned");
                    slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
                }
                drop(reader_guard);
                return Ok(true);
            }

            {
                let mut dec = slot.decompressed.lock().expect("decompressed mutex poisoned");
                for b in decompressed_batch {
                    dec.push_back(b);
                }
                if hit_eof {
                    slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
                }
            }
            drop(reader_guard);
            return Ok(true);
        }
        Ok(false)
    }
}

impl Clone for SortSpillDecompress {
    fn clone(&self) -> Self {
        Self {
            registry: Arc::clone(&self.registry),
            block_dec: SpillBlockDecompressor::new(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
        }
    }
}

impl Step for SortSpillDecompress {
    type Input = SortPhase1Event;
    type Outputs = Single<SortPhase2Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortSpillDecompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held output first.
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // 2. Pop one input event, register if SpillReady, forward all.
        if let Some(event) = ctx.input.pop() {
            let forwarded = match event {
                SortPhase1Event::SpillReady { slot, path, records_ingested_so_far } => {
                    self.registry.lock().push(RegisteredSpill { slot: Arc::clone(&slot) });
                    SortPhase2Event::SpillReady { slot, path, records_ingested_so_far }
                }
                SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far } => {
                    SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far }
                }
                SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, total_records } => {
                    SortPhase2Event::AllAnnounced { slot_count, memory_chunk_count, total_records }
                }
            };
            let _ = self.push_or_hold(ctx, forwarded);
            return Ok(StepOutcome::Progress);
        }

        // 3. Greedy slot-fill step.
        if self.try_fill_some_slot()? {
            return Ok(StepOutcome::Progress);
        }

        // 4. No fill work.
        let any_alive = self
            .snapshot_registry()
            .iter()
            .any(|slot| !slot.queue_eof.load(std::sync::atomic::Ordering::Acquire));
        if any_alive {
            return Ok(StepOutcome::Contention);
        }

        if ctx.input.is_drained() {
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests;

//! `CoalesceBytes` mid-step. `Serial + ByItemOrdinal`. Buffers
//! `DecompressedBlock` bytes from many upstream items into one combined
//! output, emitted when accumulated bytes reach a threshold. Used to
//! amortize per-item costs (notably BGZF huffman-table builds) when
//! upstream emits many small items.
//!
//! ## Why this step exists
//!
//! Legacy's `try_step_serialize` (`pipeline/bam.rs:2865-2883`)
//! accumulates *every item in a `Vec<P>` batch* into one combined
//! `SerializedBatch.data`, so by the time bytes reach the BGZF
//! compressor they're naturally in `~250 KB+` chunks — large enough
//! that each `flush()` emits several full 64 KB BGZF blocks and huffman
//! cost amortizes.
//!
//! The new pipeline's `serialize_step` is `process_ordered` — strictly
//! one-in-one-out — so each tiny position group becomes its own tiny
//! `DecompressedBlock`, which becomes its own tiny BGZF block + its own
//! huffman build. Profile diff (CODEC 8M, threads=4) showed 20× more
//! `deflate_make_huffman_code` samples in the new path than legacy.
//!
//! Inserting `CoalesceBytes` between `serialize_step` and `BgzfCompress`
//! restores the legacy property: `BgzfCompress`'s input is now a few
//! hundred KB per item, each `flush()` emits multiple full BGZF blocks,
//! and huffman cost amortizes across the block-fill.
//!
//! ## Why `Serial` (not `Parallel`)
//!
//! The accumulator state lives across `try_run` calls — one buffer
//! collecting bytes from sequential inputs. With multiple parallel
//! workers each holding their own buffer, ordinal allocation across
//! workers becomes racy and the cross-worker output ordering breaks the
//! `ByItemOrdinal` contract. `Serial` keeps the accumulator in one
//! place, mutex-protected by the framework. The work itself is just
//! `extend_from_slice` (memcpy), so single-threaded execution is fine.

use std::io;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::DecompressedBlock;

/// Default coalesce threshold: 256 KiB. Sized to fill ~4 BGZF blocks
/// per emit (BGZF max block size is ~64 KiB), matching legacy's
/// typical `SerializedBatch.data` size from `try_step_serialize`. Smaller
/// values fragment too much; larger values add latency and memory
/// without further amortizing huffman cost.
pub const DEFAULT_COALESCE_THRESHOLD_BYTES: usize = 256 * 1024;

/// Max input items consumed per `try_run` invocation. Amortizes the
/// Serial mutex acquisition.
///
/// Coalesce sees one input per upstream-emitted item. For position-group
/// workloads, those items are tiny (individual position groups, often
/// singletons), so 8 inputs per lock fills only a few KB and the
/// mutex-acquisition rate becomes the throughput cap (CPU collapsed to
/// ~1 effective core on 53M-record agilent-hs2 with a value of 8). 64
/// gives us ~64 KB worth of bytes per acquisition for typical workloads,
/// which is enough to keep the rest of the pipeline saturated. The
/// threshold check inside the loop ensures we still bound `pending`
/// memory at `threshold_bytes` regardless of batch size.
///
/// Note we do NOT mirror `GroupByPosition`'s 8 because that step's
/// per-input work (record grouping) is heavier than coalesce's
/// `extend_from_slice`; the mutex amortization budget at coalesce is
/// where we want to spend a higher count.
const MAX_BATCHES_PER_LOCK: usize = 64;

/// `Serial + ByItemOrdinal` byte coalescer for `DecompressedBlock`
/// streams. Buffers input bytes until reaching `threshold_bytes`, then
/// emits a single `DecompressedBlock` carrying the concatenated bytes
/// under a fresh monotonic ordinal.
///
/// Input items' `batch_serial` is **discarded** — coalescing decouples
/// output ordinals from input ordinals (one output corresponds to many
/// inputs). The framework's "skip `ReorderStage` on Serial producers"
/// optimization (see `core/erased.rs::TypedStep::build_output_set`)
/// elides the downstream reorder buffer because Serial pushes in
/// monotonic order by construction.
pub struct CoalesceBytes {
    threshold_bytes: usize,
    /// Self-managed monotonic output ordinal — each emitted block gets
    /// the next value. Required by the downstream `BgzfCompress`'s
    /// `BranchOrdering::ByItemOrdinal` (consumes items via
    /// `DecompressedBlock::ordinal()`).
    next_ordinal: u64,
    /// Accumulator buffer. Reset to a fresh `Vec` on each emit so the
    /// emitted block owns its bytes; the next `Vec` is sized
    /// `threshold_bytes * 2` so the typical fill cycle hits no
    /// reallocations.
    pending: Vec<u8>,
    /// Held output slot when a push was rejected by downstream
    /// backpressure.
    held: HeldSlot<Unpushed<DecompressedBlock>>,
    output_byte_limit: u64,
    name: &'static str,
}

impl CoalesceBytes {
    /// Construct a coalescer that emits when accumulated bytes reach
    /// `threshold_bytes`. Pass `output_byte_limit` for the downstream
    /// queue's byte budget (typically `BamPipelineTuning::per_step_byte_limit`).
    #[must_use]
    pub fn new(threshold_bytes: usize, output_byte_limit: u64) -> Self {
        Self {
            threshold_bytes,
            next_ordinal: 0,
            pending: Vec::with_capacity(threshold_bytes * 2),
            held: HeldSlot::new(),
            output_byte_limit,
            name: "CoalesceBytes",
        }
    }

    /// Convenience: coalescer using `DEFAULT_COALESCE_THRESHOLD_BYTES`.
    #[must_use]
    pub fn with_default_threshold(output_byte_limit: u64) -> Self {
        Self::new(DEFAULT_COALESCE_THRESHOLD_BYTES, output_byte_limit)
    }
}

impl Step for CoalesceBytes {
    type Input = DecompressedBlock;
    type Outputs = OrderedBytesSingle<DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held slot first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. If accumulator has hit the threshold, emit. Don't pull
        // more input until pending drains — keeps memory bounded at
        // ~threshold_bytes.
        if self.pending.len() >= self.threshold_bytes {
            let serial = self.next_ordinal;
            self.next_ordinal += 1;
            let bytes =
                std::mem::replace(&mut self.pending, Vec::with_capacity(self.threshold_bytes * 2));
            let out = DecompressedBlock { batch_serial: serial, bytes };
            match ctx.outputs.push(out) {
                Ok(()) => return Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        // 3. Pull up to `MAX_BATCHES_PER_LOCK` inputs per call to
        // amortize the Serial mutex acquisition. With many small
        // upstream items (each a per-position-group serialization),
        // popping one-at-a-time would make coalesce a per-byte
        // serialization point — the lock thrash dominated wall time
        // (~1 effective CPU on 4 workers). Stop early once we've
        // crossed the threshold so memory stays bounded.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(input) = ctx.input.pop() else { break };
            did_work = true;
            self.pending.extend_from_slice(&input.bytes);
            if self.pending.len() >= self.threshold_bytes {
                break;
            }
        }
        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // 4. No input this call. If upstream is drained, flush the final
        // partial block (held is empty here because step 1 early-returns
        // `Contention` on a failed retry, so control only reaches here with
        // `held` empty) and report `Finished` once nothing remains. A bounced
        // final push is parked in `held` and retried by step 1 next pass.
        if ctx.input.is_drained() {
            if !self.pending.is_empty() {
                let serial = self.next_ordinal;
                self.next_ordinal += 1;
                let bytes = std::mem::take(&mut self.pending);
                let out = DecompressedBlock { batch_serial: serial, bytes };
                if let Err(unpushed) = ctx.outputs.push(out) {
                    self.held.put(unpushed);
                }
                return Ok(StepOutcome::Progress);
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_advertises_serial_byordinal() {
        let s = CoalesceBytes::with_default_threshold(1 << 20);
        let p = s.profile();
        assert_eq!(p.name, "CoalesceBytes");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn default_threshold_is_256_kib() {
        assert_eq!(DEFAULT_COALESCE_THRESHOLD_BYTES, 256 * 1024);
    }
}

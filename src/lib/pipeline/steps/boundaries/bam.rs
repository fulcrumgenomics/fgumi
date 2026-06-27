//! `FindBamBoundaries` mid-step. `Serial + ByItemOrdinal`. Strips the BAM
//! header (on first input) and finds record boundaries within decompressed
//! BGZF block data, emitting `DecompressedBlock`s whose bytes contain only
//! complete records (4-byte `block_size` prefix + record body, repeated).
//!
//! Wraps [`super::state::BoundaryState`] — the boundary-finding state
//! machine handles header skipping, cross-block record carryover, and
//! validation. Reusing it ensures the new framework's boundary semantics
//! match the legacy pipeline's exactly.

use std::io;

use super::state::BoundaryState;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::DecompressedBlock;

/// Max inputs processed per `try_run` invocation. Amortizes the Serial
/// mutex acquisition; matches legacy `bam.rs:2050+` (`MAX_BATCHES_PER_LOCK`).
const MAX_BATCHES_PER_LOCK: usize = 8;

/// `Serial + ByItemOrdinal` boundary finder. Holds `BoundaryState` (which
/// owns the cross-block carryover buffer + header-skip flag).
pub struct FindBamBoundaries {
    state: BoundaryState,
    /// Pending output batch when we found boundaries but the push was
    /// rejected. Held across retries until pushed.
    held: HeldSlot<Unpushed<DecompressedBlock>>,
    /// Self-managed output ordinal. Incremented only when a new
    /// `DecompressedBlock` is emitted. Using the input's `batch_serial`
    /// directly would skip ordinals when the boundary state absorbs an
    /// input without producing output (header bytes, mid-record carryover);
    /// the downstream `ReorderStage` requires consecutive `0, 1, 2, …`.
    next_output_serial: u64,
    /// Set once the final-flush path has called the (non-idempotent)
    /// `state.finish()`; guards against a second call across the multi-pass
    /// completion drain.
    finalized: bool,
    output_byte_limit: u64,
}

impl FindBamBoundaries {
    /// Construct expecting the first input batch to begin with the BAM
    /// header (matches `BAM_MAGIC`). The header bytes are skipped on the
    /// first call; subsequent calls just find record boundaries.
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            state: BoundaryState::new(),
            held: HeldSlot::new(),
            next_output_serial: 0,
            finalized: false,
            output_byte_limit,
        }
    }

    /// Construct expecting the input stream to be already past the BAM
    /// header. Useful for runall-spliced sub-pipelines.
    #[must_use]
    pub fn new_no_header(output_byte_limit: u64) -> Self {
        Self {
            state: BoundaryState::new_no_header(),
            held: HeldSlot::new(),
            next_output_serial: 0,
            finalized: false,
            output_byte_limit,
        }
    }
}

impl Step for FindBamBoundaries {
    type Input = DecompressedBlock;
    type Outputs = OrderedBytesSingle<DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "FindBamBoundaries",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // Process up to `MAX_BATCHES_PER_LOCK` inputs per `try_run` to
        // amortize the Serial mutex acquisition. Each iteration pops one
        // input, runs find_boundaries, and pushes at most one output —
        // header-only or carryover-only inputs are fully absorbed and produce
        // no output (the loop `continue`s, and `did_work` still records the
        // consumed input as Progress). We stop early on push rejection (held
        // the unpushed; subsequent iterations would contend on backpressure)
        // or on input exhaustion.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(block) = ctx.input.pop() else { break };
            did_work = true;

            let boundary_batch = self.state.find_boundaries(&block.bytes)?;
            if boundary_batch.buffer.is_empty() {
                // Input fully absorbed (header/leftover); try next input.
                continue;
            }

            let serial = self.next_output_serial;
            self.next_output_serial += 1;
            let out = DecompressedBlock { batch_serial: serial, bytes: boundary_batch.buffer };
            match ctx.outputs.push(out) {
                Ok(()) => {}
                Err(unpushed) => {
                    self.held.put(unpushed);
                    // Hold off on more inputs — the held item must clear
                    // before we accept new work.
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // No input this call. If upstream is drained, flush the boundary
        // state's final partial-record buffer once (guarded by `finalized` —
        // `state.finish()` is not idempotent), emit it, and report `Finished`
        // once nothing remains. `held` is empty here (step 1 returned
        // `Contention` otherwise); a bounced final push is parked in `held`
        // for the held-drain at the top of the next pass.
        if ctx.input.is_drained() {
            if !self.finalized {
                self.finalized = true;
                if let Some(boundary_batch) = self.state.finish()? {
                    if !boundary_batch.buffer.is_empty() {
                        let serial = self.next_output_serial;
                        self.next_output_serial += 1;
                        let out = DecompressedBlock {
                            batch_serial: serial,
                            bytes: boundary_batch.buffer,
                        };
                        if let Err(unpushed) = ctx.outputs.push(out) {
                            self.held.put(unpushed);
                        }
                        return Ok(StepOutcome::Progress);
                    }
                }
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
        let s = FindBamBoundaries::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "FindBamBoundaries");
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }
}

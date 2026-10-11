//! `FrameBgzfBlocks` — cuts the sort's input slices into BGZF blocks.
//!
//! Consumes the [`ReadSlice`]s `PreadInputSlices` produces, in file order (its
//! edge is ordinal-ordered), and emits [`BgzfBlock`]s exactly as
//! `ReadBgzfBlocks` does: `batch_serial` dense from 0, EOF-marker blocks
//! skipped, `uncompressed_size` from the footer, `index: None` — except that a
//! block lying wholly inside its slice borrows it ([`RawFrame::Borrowed`])
//! instead of being copied, and only a block that straddled two slices is
//! owned. The framing
//! is [`BgzfSliceFramer`]'s, which applies the same header checks as the
//! sequential reader; a stream that ends inside a block fails with
//! `UnexpectedEof`.
//!
//! The step is reader-affine and `Serial` (not sticky: the planner holds the
//! reader worker's sticky slot). When its pop finds no slice it records that
//! in the [`InputLedger`]'s `framer_waiting`, the half of the read-streams
//! ratchet's signal that says the framer — not the device — is idle.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use fgumi_bam_io::pread::RawFrame;
use fgumi_bgzf::reader::{BgzfSliceFramer, SliceFrame, uncompressed_size};
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

use crate::pread::{InputLedger, ReadSlice};
use crate::types::BgzfBlock;

/// Counter slot index for BGZF blocks framed.
const BLOCKS: usize = 0;
/// Counter slot index for slice bytes framed.
const BYTES_READ: usize = 1;

/// `Serial` BGZF framer over ordered input slices (see the module docs).
pub struct FrameBgzfBlocks {
    framer: BgzfSliceFramer,
    ledger: Arc<InputLedger>,
    frames: Vec<SliceFrame>,
    pending: VecDeque<BgzfBlock>,
    held: HeldSlot<Unpushed<BgzfBlock>>,
    next_serial: u64,
    output_byte_limit: u64,
    finished: bool,
}

impl FrameBgzfBlocks {
    /// A framer at the start of the input, reporting to `ledger`.
    #[must_use]
    pub fn new(ledger: Arc<InputLedger>, output_byte_limit: u64) -> Self {
        Self {
            framer: BgzfSliceFramer::new(),
            ledger,
            frames: Vec::new(),
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            next_serial: 0,
            output_byte_limit,
            finished: false,
        }
    }

    /// Frame one slice into `pending`.
    fn frame(&mut self, slice: &ReadSlice) -> io::Result<()> {
        self.frames.clear();
        self.framer.push(&slice.bytes, &mut self.frames)?;
        for frame in self.frames.drain(..) {
            // A block wholly inside the slice borrows it (zero copy; the slice
            // returns to its pool when the last such block drops); only a block
            // that straddled slices was reassembled into an owned buffer.
            let bytes = match frame {
                SliceFrame::Within(r) => RawFrame::borrowed(&slice.bytes, r),
                SliceFrame::Carried(v) => RawFrame::Owned(v),
            };
            // Checked: rejects a footer claiming more than one BGZF block's
            // worth (64 KiB), so the cast below cannot truncate.
            let uncompressed_size = u32::try_from(uncompressed_size(&bytes)?)
                .expect("a checked BGZF ISIZE is at most 64 KiB");
            self.pending.push_back(BgzfBlock {
                batch_serial: self.next_serial,
                uncompressed_size,
                bytes,
                index: None,
            });
            self.next_serial += 1;
        }
        self.ledger.framed.fetch_add(slice.bytes.len() as u64, Ordering::Relaxed);
        if slice.last {
            self.framer.finish()?;
            self.finished = true;
        }
        Ok(())
    }

    /// Push the next pending block (holding it on a full queue).
    fn push_one(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        let Some(block) = self.pending.pop_front() else { return false };
        if let Err(unpushed) = ctx.outputs.push(block) {
            self.held.put(unpushed);
        }
        true
    }
}

impl Step for FrameBgzfBlocks {
    type Input = ReadSlice;
    type Outputs = OrderedBytesSingle<BgzfBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "FrameBgzfBlocks",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        Affinity::Reader
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] =
            &[CounterSpec::new("blocks", "blocks"), CounterSpec::new("bytes_read", "bytes")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        if self.push_one(ctx) {
            return Ok(StepOutcome::Progress);
        }
        if self.finished {
            return Ok(StepOutcome::Finished);
        }
        let Some(slice) = ctx.input.pop() else {
            self.ledger.framer_waiting.store(true, Ordering::Relaxed);
            if ctx.input.is_drained() {
                // The planner finished without a `last` slice only for an empty
                // input; a carried partial block is still a truncation.
                self.framer.finish()?;
                self.finished = true;
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        self.ledger.framer_waiting.store(false, Ordering::Relaxed);
        let blocks_before = self.next_serial;
        self.frame(&slice)?;
        ctx.counters.add(BLOCKS, self.next_serial - blocks_before);
        ctx.counters.add(BYTES_READ, slice.bytes.len() as u64);
        #[cfg(test)]
        AFTER_FRAME.fire(slice.bytes.len());
        drop(slice);
        self.push_one(ctx);
        Ok(StepOutcome::Progress)
    }
}

/// A test-only hook run after each slice is framed, with the slice's length.
#[cfg(test)]
pub(crate) static AFTER_FRAME: crate::test_hook::TestHook<usize> =
    crate::test_hook::TestHook::new();

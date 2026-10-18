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
//! The step is reader-affine and `Serial`, and not sticky: the planner holds
//! the reader worker's only sticky slot (`assign_sticky_owners` grants one
//! sticky step per worker). Without the sticky fast path, every dispatch is a
//! full walk of the reader worker's steps, so each dispatch pushes a batch of
//! up to `FRAME_BATCH` blocks, framing further slices as the previous one's
//! blocks drain, until the output refuses or the input runs dry. When its pop
//! finds no slice it records that in the [`InputLedger`]'s `framer_waiting`,
//! the half of the read-streams ratchet's signal that says the framer — not
//! the device — is idle.

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
    runtime::STICKY_BURST_LIMIT,
    step::{Affinity, CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

use crate::pread::{InputLedger, ReadSlice};
use crate::types::BgzfBlock;

/// Counter slot index for BGZF blocks framed.
const BLOCKS: usize = 0;
/// Counter slot index for slice bytes framed.
const BYTES_READ: usize = 1;

/// Max blocks `FrameBgzfBlocks` pushes per `try_run` dispatch.
///
/// The framer is not sticky, so it is reached only through the reader worker's
/// walk, which under the sort's drain-first scheduler restarts from the deepest
/// step on every `Progress` and so takes queued inflate work first. Pushing one
/// block per dispatch made ingest run at one block per reader-worker loop —
/// the pool starved while the framer emitted a slice's few hundred blocks one
/// walk at a time. A batch amortises the walk the same way `ReadBlocks`'
/// `ADMIT_BATCH` does on the next hop. The output edge's byte and count bounds
/// usually end the batch first; this cap bounds how long one dispatch can keep
/// the reader worker when the consumers drain the edge as fast as it fills.
///
/// It is the runtime's sticky burst bound ([`STICKY_BURST_LIMIT`]), the most
/// consecutive calls a sticky step such as `ReadBgzfBlocks` gets before its
/// worker walks its other steps; a one-block-per-call sticky reader therefore
/// emits at most this many blocks per visit, and so does this batch. The value
/// is the runtime's, not a copy, so the two cannot drift apart.
pub(crate) const FRAME_BATCH: usize = STICKY_BURST_LIMIT;

/// `Serial` BGZF framer over ordered input slices (see the module docs).
pub struct FrameBgzfBlocks {
    framer: BgzfSliceFramer,
    ledger: Arc<InputLedger>,
    frames: Vec<SliceFrame>,
    pending: VecDeque<BgzfBlock>,
    held: HeldSlot<Unpushed<BgzfBlock>>,
    next_serial: u64,
    output_byte_limit: u64,
    /// Max blocks offered per dispatch: `FRAME_BATCH`, except in tests that
    /// need a cap below the output edge's slot count to observe it.
    batch: usize,
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
            batch: FRAME_BATCH,
            finished: false,
        }
    }

    /// This framer with a per-dispatch cap of `batch` blocks instead of
    /// `FRAME_BATCH`. `FRAME_BATCH` equals the output edge's slot count, so a
    /// hand-driven probe fills the edge before it reaches the cap and cannot
    /// tell a capped batch from an uncapped one; a smaller cap can.
    #[cfg(test)]
    pub(crate) fn with_batch(mut self, batch: usize) -> Self {
        assert!(batch > 0, "a zero cap would never push");
        self.batch = batch;
        self
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

    /// Push pending blocks until `pending` is empty, the output refuses one
    /// (it is held for the next dispatch), or `budget` blocks have been
    /// offered. Returns how many blocks were offered, and whether the output
    /// refused.
    fn push_pending(&mut self, ctx: &mut StepCtx<'_, Self>, budget: usize) -> (usize, bool) {
        let mut offered = 0;
        while offered < budget {
            let Some(block) = self.pending.pop_front() else { break };
            offered += 1;
            if let Err(unpushed) = ctx.outputs.push(block) {
                self.held.put(unpushed);
                return (offered, true);
            }
        }
        (offered, false)
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
        // Batched (see `FRAME_BATCH`): push what is pending, then frame the
        // next slice and keep pushing, until the output refuses, the batch is
        // spent, or no slice is ready. A popped slice is progress even when it
        // framed no whole block (it freed a slot on the input edge).
        let mut offered = 0;
        let mut popped_any = false;
        loop {
            let (n, refused) = self.push_pending(ctx, self.batch - offered);
            offered += n;
            if refused || offered >= self.batch {
                return Ok(StepOutcome::Progress);
            }
            let progressed = offered > 0 || popped_any;
            if self.finished {
                return Ok(if progressed { StepOutcome::Progress } else { StepOutcome::Finished });
            }
            let Some(slice) = ctx.input.pop() else {
                self.ledger.framer_waiting.store(true, Ordering::Relaxed);
                if progressed {
                    return Ok(StepOutcome::Progress);
                }
                if ctx.input.is_drained() {
                    // The planner finished without a `last` slice only for an
                    // empty input; a carried partial block is still a truncation.
                    self.framer.finish()?;
                    self.finished = true;
                    return Ok(StepOutcome::Finished);
                }
                return Ok(StepOutcome::NoProgress);
            };
            popped_any = true;
            self.ledger.framer_waiting.store(false, Ordering::Relaxed);
            let blocks_before = self.next_serial;
            self.frame(&slice)?;
            ctx.counters.add(BLOCKS, self.next_serial - blocks_before);
            ctx.counters.add(BYTES_READ, slice.bytes.len() as u64);
            #[cfg(test)]
            AFTER_FRAME.fire(slice.bytes.len());
            drop(slice);
        }
    }
}

/// A test-only hook run after each slice is framed, with the slice's length.
#[cfg(test)]
pub(crate) static AFTER_FRAME: crate::test_hook::TestHook<usize> =
    crate::test_hook::TestHook::new();

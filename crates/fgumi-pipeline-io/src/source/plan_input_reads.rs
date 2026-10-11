//! `PlanInputReads` — the planner of the sort's chain-native input reads.
//!
//! A `Serial`, sticky, reader-affine source that decides which byte ranges of
//! the input BAM to read next and emits them as [`ReadRequest`]s for
//! `PreadInputSlices` (a [`crate::pread::PreadSlices`]) to perform; `FrameBgzfBlocks`
//! then cuts the ordered slices into BGZF blocks. Each fill is
//! [`FILL_BYTES`], split into [`ReadStreamsPolicy::slices_for`] equal slices.
//!
//! Read-ahead is bounded by its own ledger, not by the queue limits: the
//! planner issues a fill only while `issued − framed` leaves room for one more
//! fill within the lookahead ([`LOOKAHEAD_BYTES`], two fills), so at most
//! two fills are outstanding or waiting to be framed.
//!
//! The planner also drives the `--read-streams auto` ratchet: on every visit
//! it reports whether the device is the bottleneck
//! ([`fgumi_bam_io::pread::starved_predicate`] over the shared
//! [`InputLedger`]) — every permitted byte outstanding *and* the framer
//! waiting for its head slice. A slow consumer never ratchets, because the
//! framer is then busy rather than waiting.

use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use fgumi_bam_io::pread::{
    FILL_BYTES, PositionalSource, ReadStreamsPolicy, slice_ranges, starved_predicate,
};
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

use crate::pread::{InputLedger, ReadRequest, ReadTarget};

/// Fills the planner may run ahead of the framer beyond the one being framed.
pub const PHASE1_LOOKAHEAD_FILLS: u64 = 1;

/// The input lookahead: `(PHASE1_LOOKAHEAD_FILLS + 1) × FILL_BYTES` = 8 MiB.
pub const LOOKAHEAD_BYTES: u64 = (PHASE1_LOOKAHEAD_FILLS + 1) * FILL_BYTES as u64;

/// Counter slot index for fills planned.
const FILLS: usize = 0;
/// Counter slot index for read requests emitted.
const REQUESTS: usize = 1;

/// `Serial` planner of the input BAM's positional reads (see the module docs).
pub struct PlanInputReads {
    source: Arc<dyn PositionalSource>,
    len: u64,
    next_offset: u64,
    next_seq: u32,
    next_ordinal: u64,
    policy: Arc<ReadStreamsPolicy>,
    ledger: Arc<InputLedger>,
    eligible_clones: usize,
    fills_since_observe: u32,
    pending: VecDeque<ReadRequest>,
    held: HeldSlot<Unpushed<ReadRequest>>,
    output_byte_limit: u64,
}

impl PlanInputReads {
    /// Plan reads of the file at `path`. `eligible_clones` is the number of
    /// `PreadInputSlices` clones that can read at once (more slices per fill
    /// than that only queue).
    ///
    /// # Errors
    /// I/O errors opening the file or reading its length.
    pub fn open(
        path: &Path,
        policy: Arc<ReadStreamsPolicy>,
        ledger: Arc<InputLedger>,
        eligible_clones: usize,
        output_byte_limit: u64,
    ) -> io::Result<Self> {
        let file = File::open(path)?;
        let len = PositionalSource::byte_len(&file)?;
        Ok(Self::from_source(
            Arc::new(file),
            len,
            policy,
            ledger,
            eligible_clones,
            output_byte_limit,
        ))
    }

    /// Plan reads of `len` bytes of `source`.
    #[must_use]
    pub fn from_source(
        source: Arc<dyn PositionalSource>,
        len: u64,
        policy: Arc<ReadStreamsPolicy>,
        ledger: Arc<InputLedger>,
        eligible_clones: usize,
        output_byte_limit: u64,
    ) -> Self {
        Self {
            source,
            len,
            next_offset: 0,
            next_seq: 0,
            next_ordinal: 0,
            policy,
            ledger,
            eligible_clones,
            fills_since_observe: 0,
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            output_byte_limit,
        }
    }

    /// Queue one fill's requests and charge it to the ledger.
    fn issue_fill(&mut self) {
        let remaining = self.len - self.next_offset;
        let want = usize::try_from(remaining).map_or(FILL_BYTES, |r| r.min(FILL_BYTES));
        let slices = self.policy.slices_for(want, self.eligible_clones);
        for (offset, len) in slice_ranges(self.next_offset, want, slices) {
            self.pending.push_back(ReadRequest {
                ordinal: self.next_ordinal,
                stream: 0,
                seq: self.next_seq,
                source: Arc::clone(&self.source),
                offset,
                len,
                last: offset + u64::from(len) == self.len,
                target: ReadTarget::Input,
            });
            self.next_ordinal += 1;
            self.next_seq += 1;
        }
        self.next_offset += want as u64;
        self.ledger.issued.fetch_add(want as u64, Ordering::Relaxed);
        self.fills_since_observe += 1;
    }
}

impl Step for PlanInputReads {
    type Input = ();
    type Outputs = OrderedBytesSingle<ReadRequest>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "PlanInputReads",
            kind: StepKind::Serial,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        Affinity::Reader
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] =
            &[CounterSpec::new("fills", "fills"), CounterSpec::new("requests", "requests")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        if let Some(req) = self.pending.pop_front() {
            if let Err(unpushed) = ctx.outputs.push(req) {
                self.held.put(unpushed);
            }
            ctx.counters.add(REQUESTS, 1);
            return Ok(StepOutcome::Progress);
        }

        let issued = self.ledger.issued.load(Ordering::Relaxed);
        let landed = self.ledger.landed.load(Ordering::Relaxed);
        let framed = self.ledger.framed.load(Ordering::Relaxed);
        let waiting = self.ledger.framer_waiting.load(Ordering::Relaxed);
        let starved = starved_predicate(issued, landed, LOOKAHEAD_BYTES, waiting);
        self.policy.observe(starved, std::mem::take(&mut self.fills_since_observe));

        if self.next_offset == self.len {
            return Ok(StepOutcome::Finished);
        }
        if issued.saturating_sub(framed) + FILL_BYTES as u64 > LOOKAHEAD_BYTES {
            return Ok(StepOutcome::NoProgress);
        }
        self.issue_fill();
        ctx.counters.add(FILLS, 1);
        #[cfg(any(test, feature = "test-utils"))]
        test_hooks::after_fill(self.next_offset == self.len);
        Ok(StepOutcome::Progress)
    }
}

/// A hook run after each planned fill, on the planner's thread, with whether
/// the fill was the last (test support: lets a test observe the process at
/// read time without timing).
#[cfg(any(test, feature = "test-utils"))]
pub mod test_hooks {
    use crate::test_hook::TestHook;

    static AFTER_FILL: TestHook<bool> = TestHook::new();

    /// Install `f`, called with `last` after each fill (one per process: nextest
    /// runs each test in its own).
    pub fn set(f: impl Fn(bool) + Send + Sync + 'static) {
        AFTER_FILL.set(f);
    }

    /// Remove the hook.
    pub fn clear() {
        AFTER_FILL.clear();
    }

    pub(super) fn after_fill(last: bool) {
        AFTER_FILL.fire(last);
    }
}

//! `SortSpillDecompress` — ingests the merge's spill read slices and
//! decompresses their frames into the slots' bounded queues.
//!
//! A `Parallel` leaf (`Outputs = ()`): its input is the [`ReadSlice`]s
//! `PreadSpillSlices` read for `SpillReadPlanner`; its output is each slot's
//! decompressed FIFO, a side channel the merge pops. One `try_run` does one
//! of three things:
//!
//! 1. **Ingest** one slice: the target slot reorders it by its per-slot slice
//!    sequence, parses every frame of each in-order slice into its raw stash
//!    ([`SortMergeSlot::bp_ingest_slice`]), and this clone then serves that
//!    slot once.
//! 2. **Serve** one block: claim a stash head ([`SortMergeSlot::bp_claim_raw`]),
//!    decompress it, publish it ([`SortMergeSlot::bp_insert_drain_finalize`]),
//!    and wake the merge if it awaits the slot. Slots are scanned demand first
//!    — the merge's awaited slot, then a rotation from a cursor the clones share —
//!    using lock-free mirrors only, so a scan never takes a slot lock.
//! 3. **Drain**: move a block waiting in a slot's reorder buffer into its FIFO
//!    once the consumer has made room, and finalize `queue_eof` — the path that
//!    delivers a last block inserted while the FIFO was full.
//!
//! # Ownership of a raw block, and why the merge cannot be left behind
//!
//! A raw block is owned by its **slot** from its parse
//! ([`SortMergeSlot::bp_commit_read`] reserves it as in flight) until its
//! decompressed bytes are inserted; any clone — or the merge itself — may
//! claim it in between. Liveness:
//!
//! - Someone keeps visiting: a clone reports `Finished` only once its input is
//!   drained *and* every registered slot is `queue_eof`; a stashed block holds
//!   `in_flight > 0`, which blocks `queue_eof`, so no clone leaves while one
//!   exists. Each `Progress` of this zero-output `Parallel` step cascades a
//!   wake to the pool, so a peer can claim what an ingest exposed.
//! - The awaited slot is always claimable: the merge stalls only on an empty
//!   FIFO, so the FIFO cap admits its head, and the head that is the reorder
//!   front is admitted unconditionally (the escape in `bp_claim_raw`).
//! - Backpressure composes without a queue cycle: claims stop at the FIFO cap,
//!   the merge's pops make room, reads stop at the planner's allowance, and
//!   claims drain the allowance.
//!
//! # Admission
//!
//! Work is checked before the phase-2 permit is taken: a clone with no slice
//! to ingest and no claimable or drainable slot (read from the slots'
//! lock-free mirrors) takes no permit and reports `NoProgress` (or `Finished`).

use std::collections::HashMap;
use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use fgumi_sort::{SortMergeSlot, SpillBlockDecompressor};
use parking_lot::RwLock;

use super::supply_ledger::{SpillSupply, SupplyLedger};
use crate::pread::{ReadSlice, ReadTarget};
use fgumi_pipeline_core::{
    PhaseCap,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Default reorder-window byte budget substituted when the caller's
/// `output_byte_limit` is `0`. A claim reads `window_budget == 0` as *no
/// bound*, so a resolved-to-zero budget is normalized to this cap rather than
/// removing the only byte cap on decompressed stragglers.
const DEFAULT_REORDER_WINDOW_BYTES: u64 = 256 * 1024 * 1024;

/// The slots seen so far (append-only), and their index by file id.
#[derive(Default)]
struct Registry {
    slots: Vec<Arc<SortMergeSlot>>,
    by_file_id: HashMap<u32, usize>,
}

/// `Parallel` spill-slice ingest + decompress step (see the module docs).
pub struct SortSpillDecompress {
    registry: Arc<RwLock<Registry>>,
    block_dec: SpillBlockDecompressor,
    /// Phase-2 admission cap shared with the other merge-phase pool steps
    /// (`None` = uncapped).
    cap: Option<Arc<PhaseCap>>,
    /// Per-slot reorder-window byte budget for claims.
    window_budget: u64,
    /// The sort's merge-wide demand: deliveries to the awaited slot wake the
    /// merge; the merge's hot set is scanned first.
    demand: Arc<fgumi_sort::MergeDemand>,
    /// The supply ledger shared with the read planner: landing a slice
    /// releases the planner's outstanding-read bound.
    ledger: Arc<SupplyLedger>,
    /// The rotation start, shared by every clone: each scan takes the next
    /// one, so concurrent clones start at different slots.
    cursor: Arc<std::sync::atomic::AtomicUsize>,
}

impl SortSpillDecompress {
    /// A step whose claims are bounded by a per-slot reorder window of
    /// `output_byte_limit` bytes, sharing `supply` (merge demand and ledger)
    /// with the read planner and the merge.
    #[must_use]
    pub fn new(output_byte_limit: u64, supply: &SpillSupply) -> Self {
        let window_budget =
            if output_byte_limit == 0 { DEFAULT_REORDER_WINDOW_BYTES } else { output_byte_limit };
        Self {
            registry: Arc::new(RwLock::new(Registry::default())),
            block_dec: SpillBlockDecompressor::new(),
            cap: None,
            window_budget,
            demand: Arc::clone(&supply.demand),
            ledger: Arc::clone(&supply.ledger),
            cursor: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
        }
    }

    /// Share the phase-2 admission cap with the other merge-phase pool steps
    /// (`--merge-threads`; `None` = uncapped).
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    /// The merge demand this step notifies (tests: a clone must carry it).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn merge_demand_for_test(&self) -> &Arc<fgumi_sort::MergeDemand> {
        &self.demand
    }

    /// Register `slot` (test support).
    #[cfg(test)]
    pub(crate) fn register_for_test(&self, slot: &Arc<SortMergeSlot>) {
        self.register(slot);
    }

    fn register(&self, slot: &Arc<SortMergeSlot>) {
        if self.registry.read().by_file_id.contains_key(&slot.file_id) {
            return;
        }
        let mut r = self.registry.write();
        if !r.by_file_id.contains_key(&slot.file_id) {
            let i = r.slots.len();
            r.by_file_id.insert(slot.file_id, i);
            r.slots.push(Arc::clone(slot));
        }
    }

    /// Lock-free: `slot` has a stash head a claim can admit, or a block a
    /// drain can move, or an EOF a drain can finalize.
    fn has_work(slot: &SortMergeSlot) -> bool {
        if slot.queue_eof() {
            return false;
        }
        let stashed = slot.stash_len_relaxed();
        let in_flight = slot.in_flight();
        let room = slot.fifo_len_relaxed() < slot.fifo_cap() as usize;
        let reorder_len = slot.reorder_len_relaxed();
        // Without FIFO room only the front escape admits a claim, and the stash
        // head is the reorder front only when no claimed block is outstanding
        // (`in_flight` counts exactly the stash then) and none waits in reorder.
        let front = in_flight == stashed && reorder_len == 0;
        let claimable = stashed > 0 && (room || front);
        let drainable = room && reorder_len > 0;
        let finalizable = slot.reader_eof() && in_flight == 0 && reorder_len == 0;
        claimable || drainable || finalizable
    }

    /// The registry indices this clone visits, in order, and whether any of
    /// them has work: demand first (the merge's hot set,
    /// [`fgumi_sort::MergeDemand::hot_ids`], unless EOF), then a rotation from
    /// the shared cursor over the slots with work. One pass over the slots'
    /// atomics, never a slot lock; returns indices so a scan clones no slot
    /// handle it does not serve.
    pub(crate) fn scan_order(&self) -> (Vec<usize>, bool) {
        let r = self.registry.read();
        let k = r.slots.len();
        let mut order = Vec::new();
        let mut any_work = false;
        for id in self.demand.hot_ids().as_slice() {
            if let Some(&i) = r.by_file_id.get(id)
                && !order.contains(&i)
                && !r.slots[i].queue_eof()
            {
                any_work |= Self::has_work(&r.slots[i]);
                order.push(i);
            }
        }
        let demand = order.len();
        if k > 0 {
            let start = self.cursor.fetch_add(1, Ordering::Relaxed) % k;
            for j in 0..k {
                let i = (start + j) % k;
                if !order[..demand].contains(&i) && Self::has_work(&r.slots[i]) {
                    any_work = true;
                    order.push(i);
                }
            }
        }
        (order, any_work)
    }

    /// The registered slot at index `i`.
    fn slot_at(&self, i: usize) -> Arc<SortMergeSlot> {
        Arc::clone(&self.registry.read().slots[i])
    }

    fn all_eof(&self) -> bool {
        self.registry.read().slots.iter().all(|s| s.queue_eof())
    }

    /// Claim, decompress and publish one block of `slot`
    /// ([`fgumi_sort::ClaimedBlock::decompress_and_publish`]). Returns whether
    /// a block was claimed.
    fn serve_one(&mut self, slot: &SortMergeSlot) -> io::Result<bool> {
        let Some(block) = slot.bp_claim_raw(self.window_budget) else { return Ok(false) };
        self.ledger.sub_stash(block.frame.len() as u64, false);
        match block.decompress_and_publish(&mut self.block_dec) {
            Ok(progressed) => {
                if progressed {
                    self.notify(slot);
                }
                Ok(true)
            }
            Err(e) => {
                self.notify(slot);
                Err(e)
            }
        }
    }

    /// Ingest one slice and serve its slot once.
    fn ingest(&mut self, slice: ReadSlice) -> io::Result<()> {
        let ReadTarget::Spill { slot, class } = slice.target else {
            return Err(io::Error::other("SortSpillDecompress received an input-path slice"));
        };
        self.register(&slot);
        let len = slice.bytes.len() as u64;
        let out = match slot.bp_ingest_slice(slice.seq, slice.bytes, slice.last) {
            Ok(out) => out,
            Err(e) => {
                self.notify(&slot);
                return Err(e);
            }
        };
        self.ledger.land(class, len);
        self.ledger.add_stash(out.stashed_bytes, out.frames as u64);
        if out.hit_eof && slot.queue_eof() {
            self.notify(&slot);
        }
        self.serve_one(&slot)?;
        Ok(())
    }

    fn run_once(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // A slice to ingest is work; otherwise one scan both decides whether
        // there is work (before a permit is taken) and orders the serve.
        let order = if ctx.input.is_empty() {
            let (order, any_work) = self.scan_order();
            if !any_work {
                if ctx.input.is_drained() && self.all_eof() {
                    return Ok(StepOutcome::Finished);
                }
                ctx.input.note_empty_poll();
                return Ok(StepOutcome::NoProgress);
            }
            Some(order)
        } else {
            None
        };
        let cap = self.cap.take();
        let outcome = self.run_admitted(ctx, cap.as_deref(), order);
        self.cap = cap;
        outcome
    }

    fn run_admitted(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
        order: Option<Vec<usize>>,
    ) -> io::Result<StepOutcome> {
        let _permit = match cap {
            Some(cap) => match cap.try_acquire() {
                Some(p) => Some(p),
                None => return Ok(StepOutcome::Capped),
            },
            None => None,
        };
        if let Some(slice) = ctx.input.pop() {
            self.ingest(slice)?;
            return Ok(StepOutcome::Progress);
        }
        let order = order.unwrap_or_else(|| self.scan_order().0);
        for i in order {
            let slot = self.slot_at(i);
            if self.serve_one(&slot)? {
                return Ok(StepOutcome::Progress);
            }
            if slot.bp_drain_and_finalize() {
                self.notify(&slot);
                return Ok(StepOutcome::Progress);
            }
        }
        Ok(StepOutcome::Contention)
    }

    /// Tell the merge that `slot` received a block, EOF or failure (after the
    /// slot's `decompressed` mutex is released).
    #[inline]
    fn notify(&self, slot: &SortMergeSlot) {
        self.demand.notify_delivered(slot.file_id);
    }
}

impl Clone for SortSpillDecompress {
    fn clone(&self) -> Self {
        Self {
            registry: Arc::clone(&self.registry),
            block_dec: SpillBlockDecompressor::new(),
            cap: self.cap.clone(),
            window_budget: self.window_budget,
            demand: self.demand.clone(),
            ledger: self.ledger.clone(),
            cursor: Arc::clone(&self.cursor),
        }
    }
}

impl Step for SortSpillDecompress {
    type Input = ReadSlice;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortSpillDecompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        self.run_once(ctx)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn phase_cap(&self) -> Option<&PhaseCap> {
        self.cap.as_deref()
    }
}

#[cfg(test)]
mod tests;

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
//!    — the merge's awaited, predicted and frontier slots, then a rotation from
//!    a cursor the clones share —
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
//! - The merge serves itself: a stalled merge claims and decompresses its
//!   awaited slot's stashed blocks on its own thread before it parks, so a
//!   block that has been read never waits for a worker to be scheduled.
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
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use fgumi_sort::{SortMergeSlot, SpillBlockDecompressor};
use parking_lot::RwLock;

use super::supply_ledger::{SpillSupply, SupplyLedger};
use crate::pread::{ReadSlice, ReadTarget};
use fgumi_pipeline_core::{
    PhaseCap, StepCap, StepCapRef,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Default reorder-window byte budget substituted when the caller's
/// `output_byte_limit` is `0`. A claim reads `window_budget == 0` as *no
/// bound*, so a resolved-to-zero budget is normalized to this cap rather than
/// removing the only byte cap on decompressed stragglers.
const DEFAULT_REORDER_WINDOW_BYTES: u64 = 256 * 1024 * 1024;

/// The per-slot reorder-window byte budget for claims from an
/// `output_byte_limit` (zero is normalized to [`DEFAULT_REORDER_WINDOW_BYTES`]).
/// Shared by the decompress step and the merge's self-serve, which claim from
/// the same stashes under the same window.
pub(crate) fn reorder_window_budget(output_byte_limit: u64) -> u64 {
    if output_byte_limit == 0 { DEFAULT_REORDER_WINDOW_BYTES } else { output_byte_limit }
}

/// The live slots, and their index by slot id. Retired slots are removed
/// (`drop_retired` compacts the vector and bumps `generation`), so the registry
/// stays as large as the live slot set. A scan returns indices together with
/// the generation they belong to; an index from an older generation is never
/// resolved (`slot_at`), so compaction by one clone cannot redirect another
/// clone's pending scan to a different slot.
#[derive(Default)]
struct Registry {
    slots: Vec<Arc<SortMergeSlot>>,
    by_slot_id: HashMap<u32, usize>,
    generation: u64,
}

/// `Parallel` spill-slice ingest + decompress step (see the module docs).
pub struct SortSpillDecompress {
    registry: Arc<RwLock<Registry>>,
    block_dec: SpillBlockDecompressor,
    /// Phase-2 admission cap shared with the other merge-phase pool steps
    /// (`None` = uncapped).
    cap: Option<StepCap>,
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
    /// The consolidation demand, and the flag the read planner sets at
    /// `AllAnnounced`: deliveries also notify the consolidation merge, and its
    /// hot set joins the scan, only before that (`None`: no consolidation tail).
    consolidation: Option<(Arc<fgumi_sort::MergeDemand>, Arc<AtomicBool>)>,
    /// The planner's retire epoch (shared through [`SpillSupply`]), and the
    /// value this clone last acted on.
    retire_epoch: Arc<AtomicU64>,
    seen_epoch: u64,
    /// This clone's scan buffer, reused across polls.
    scan_scratch: Vec<usize>,
}

impl SortSpillDecompress {
    /// A step whose claims are bounded by a per-slot reorder window of
    /// `output_byte_limit` bytes, sharing `supply` (merge demand and ledger)
    /// with the read planner and the merge.
    #[must_use]
    pub fn new(output_byte_limit: u64, supply: &SpillSupply) -> Self {
        let window_budget = reorder_window_budget(output_byte_limit);
        Self {
            registry: Arc::new(RwLock::new(Registry::default())),
            block_dec: SpillBlockDecompressor::new(),
            cap: None,
            window_budget,
            demand: Arc::clone(&supply.demand),
            ledger: Arc::clone(&supply.ledger),
            cursor: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            consolidation: None,
            retire_epoch: Arc::clone(&supply.retire_epoch),
            seen_epoch: 0,
            scan_scratch: Vec::new(),
        }
    }

    /// Also notify `demand` (the consolidation merge's) of deliveries, and scan
    /// its hot set first, until `all_announced` is set (the read planner sets it
    /// at `AllAnnounced`; no consolidation runs after it, so phase-2 deliveries
    /// pay nothing for it).
    #[must_use]
    pub fn with_consolidation_demand(
        mut self,
        demand: Arc<fgumi_sort::MergeDemand>,
        all_announced: Arc<AtomicBool>,
    ) -> Self {
        self.consolidation = Some((demand, all_announced));
        self
    }

    /// The live consolidation demand: set, and `AllAnnounced` not yet seen.
    fn consolidation_demand(&self) -> Option<&fgumi_sort::MergeDemand> {
        self.consolidation
            .as_ref()
            .filter(|(_, announced)| !announced.load(Ordering::Acquire))
            .map(|(d, _)| &**d)
    }

    /// Drop retired slots from the registry when the planner's retire epoch
    /// moved since this clone last looked (one `Acquire` load otherwise), so
    /// their files close.
    fn drop_retired(&mut self) {
        let epoch = self.retire_epoch.load(Ordering::Acquire);
        if epoch == self.seen_epoch {
            return;
        }
        self.seen_epoch = epoch;
        let mut r = self.registry.write();
        let Registry { slots, by_slot_id, generation } = &mut *r;
        let before = slots.len();
        slots.retain(|slot| !slot.is_retired());
        if slots.len() != before {
            by_slot_id.clear();
            by_slot_id.extend(slots.iter().enumerate().map(|(i, s)| (s.slot_id(), i)));
            *generation += 1;
        }
    }

    /// The live registered slots' ids (tests).
    #[cfg(test)]
    pub(crate) fn registered_slot_ids_for_test(&self) -> Vec<u32> {
        self.registry.read().slots.iter().map(|s| s.slot_id()).collect()
    }

    /// Admit through `cap` (`None` = uncapped): the phase-2 cap shared with
    /// the other merge-phase pool steps (`--merge-threads`), or a phased one
    /// that switches from phase 1 to phase 2 at `AllAnnounced`.
    #[must_use]
    pub fn with_step_cap(mut self, cap: impl Into<Option<StepCap>>) -> Self {
        self.cap = cap.into();
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
        let id = slot.slot_id();
        if self.registry.read().by_slot_id.contains_key(&id) || slot.is_retired() {
            return;
        }
        let mut r = self.registry.write();
        if !r.by_slot_id.contains_key(&id) {
            let i = r.slots.len();
            r.by_slot_id.insert(id, i);
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

    /// Fill `order` with the registry indices this clone visits, in order, and
    /// return the registry generation they belong to and whether any of them
    /// has work: demand first (the merge's hot set,
    /// [`fgumi_sort::MergeDemand::hot_ids`], unless EOF), then a rotation from
    /// the shared cursor over the slots with work. One pass over the slots'
    /// atomics, never a slot lock; returns indices, into the caller's reused
    /// buffer, so a scan neither allocates nor clones a slot handle it does
    /// not serve.
    pub(crate) fn scan_order(&self, order: &mut Vec<usize>) -> (u64, bool) {
        order.clear();
        let r = self.registry.read();
        let k = r.slots.len();
        let mut any_work = false;
        let merge_hot = self.demand.hot_ids();
        let cons_hot = self.consolidation_demand().map(fgumi_sort::MergeDemand::hot_ids);
        let hot = merge_hot
            .as_slice()
            .iter()
            .chain(cons_hot.iter().flat_map(fgumi_sort::HotIds::as_slice));
        for id in hot {
            if let Some(&i) = r.by_slot_id.get(id)
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
        (r.generation, any_work)
    }

    /// The registered slot at index `i` of registry `generation`; `None` once
    /// a compaction moved on from that generation (the caller rescans next
    /// poll).
    fn slot_at(&self, generation: u64, i: usize) -> Option<Arc<SortMergeSlot>> {
        let r = self.registry.read();
        (r.generation == generation).then(|| Arc::clone(&r.slots[i]))
    }

    fn all_eof(&self) -> bool {
        self.registry.read().slots.iter().all(|s| s.queue_eof())
    }

    /// Claim, decompress and publish one block of `slot`
    /// ([`fgumi_sort::ClaimedBlock::decompress_and_publish`]). Returns whether
    /// a block was claimed.
    fn serve_one(&mut self, slot: &SortMergeSlot) -> io::Result<bool> {
        let Some(block) = slot.bp_claim_raw(self.window_budget) else { return Ok(false) };
        self.ledger.sub_stash(block.charge, false);
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
        self.drop_retired();
        // A slice to ingest is work; otherwise one scan both decides whether
        // there is work (before a permit is taken) and orders the serve.
        let mut order = std::mem::take(&mut self.scan_scratch);
        let scanned = if ctx.input.is_empty() {
            let (generation, any_work) = self.scan_order(&mut order);
            if !any_work {
                self.scan_scratch = order;
                if ctx.input.is_drained() && self.all_eof() {
                    return Ok(StepOutcome::Finished);
                }
                ctx.input.note_empty_poll();
                return Ok(StepOutcome::NoProgress);
            }
            Some(generation)
        } else {
            None
        };
        let cap = self.cap.take();
        let outcome =
            self.run_admitted(ctx, cap.as_ref().and_then(StepCap::active), scanned, &mut order);
        self.cap = cap;
        self.scan_scratch = order;
        outcome
    }

    fn run_admitted(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
        scanned: Option<u64>,
        order: &mut Vec<usize>,
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
        let generation = match scanned {
            Some(generation) => generation,
            None => self.scan_order(order).0,
        };
        for &i in order.iter() {
            let Some(slot) = self.slot_at(generation, i) else { break };
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
        self.demand.notify_delivered(slot.slot_id());
        if let Some(d) = self.consolidation_demand() {
            d.notify_delivered(slot.slot_id());
        }
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
            consolidation: self.consolidation.clone(),
            retire_epoch: Arc::clone(&self.retire_epoch),
            seen_epoch: self.seen_epoch,
            scan_scratch: Vec::new(),
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

    fn step_cap(&self) -> Option<StepCapRef<'_>> {
        self.cap.as_ref().map(StepCap::view)
    }
}

#[cfg(test)]
mod tests;

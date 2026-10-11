//! `SpillReadPlanner` — decides which spill bytes the merge reads next.
//!
//! A pool `Serial` step between the spill writer and the merge. It owns the
//! table of merge slots, forwards the phase events to `SortMerge` on its second
//! output, and on its first emits byte-range [`ReadRequest`]s that
//! `PreadSpillSlices` performs; `SortSpillDecompress` parses the landed slices
//! into each slot's raw stash and decompresses them.
//!
//! Read-ahead is bounded by its own byte accounting, not by queue limits,
//! against a `ReadAheadBudget` resolved at `AllAnnounced` from the actual
//! slot count. Each slot is charged what it holds (requested bytes plus
//! `SortMergeSlot::stash_bytes`). Every slot is kept topped up to its cold
//! allowance in cold fills (cold slots are visited in a rotation); the slots
//! the merge needs next (the hot set) top up further, in 4 MiB fills toward
//! the hot allowance, only from the budget's shared pool, which every byte a
//! slot holds above its cold allowance — hot or not, and idle slice buffers —
//! is charged to. So a slot that leaves the hot set keeps its charge until the
//! merge drains it, and read-ahead never exceeds `k × cold + pool`. Outstanding
//! slices are capped at twice the phase-2 thread count (a read occupies a pool
//! worker), and outstanding cold bytes at `COLD_INFLIGHT_BYTES`, so hot reads
//! never queue behind a wall of cold ones; a hot slot's reads skip that cap, so
//! the awaited slot can always read within its cold terms. The planner reads
//! only lock-free state ([`fgumi_sort::MergeDemand`]'s published slots, each
//! slot's mirrors and the slice pool's gauge).
//!
//! The planner retires once every spill byte is requested, so the FIFO caps
//! are frozen for the merge's tail (see `try_run`).
//!
//! Slices per fill adopt the input's `--read-streams` count
//! ([`ReadStreamsPolicy::slices_for`]); the planner never moves the ratchet.

use std::collections::{HashMap, VecDeque};
use std::io;
use std::sync::Arc;

use fgumi_bam_io::pread::{PositionalSource, ReadStreamsPolicy, SliceBufferPool, slice_ranges};
use fgumi_pipeline_core::{
    HeldSlot, Unpushed,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};
use fgumi_sort::{MergeDemand, SortMergeSlot};

use super::protocol::{SortPhase1Event, SortPhase2Event};
use super::read_ahead_budget::{
    COLD_INFLIGHT_BYTES, FIFO_CAP_COLD, FIFO_CAP_HOT, ReadAheadBudget, fifo_charge,
};
use super::supply_ledger::{SpillSupply, SupplyLedger};
use crate::pread::{ReadRequest, ReadTarget, SpillClass};

/// Counter slot: fills issued.
const FILLS: usize = 0;
/// Counter slot: bytes issued.
const BYTES_ISSUED: usize = 1;

/// One registered slot and the planner's cursor into it.
struct PlannedSlot {
    slot: Arc<SortMergeSlot>,
    source: Arc<dyn PositionalSource>,
    next_offset: u64,
    next_seq: u32,
    eof_issued: bool,
}

/// `Serial` planner of the merge's spill reads (see the module docs).
pub struct SpillReadPlanner {
    slots: Vec<PlannedSlot>,
    by_file_id: HashMap<u32, usize>,
    demand: Arc<MergeDemand>,
    policy: Arc<ReadStreamsPolicy>,
    ledger: Arc<SupplyLedger>,
    /// The spill slice pool: its idle buffers are charged to the budget's pool.
    slices: Arc<SliceBufferPool>,
    total_memory: u64,
    max_inflight_slices: u32,
    cold_inflight_bytes: u64,
    eligible_clones: usize,
    budget: Option<ReadAheadBudget>,
    /// The hot set of the previous pass, demoted to the cold FIFO cap when
    /// it leaves the set.
    hot_prev: Vec<usize>,
    cursor: usize,
    next_ordinal: u64,
    outbox: VecDeque<ReadRequest>,
    held_read: HeldSlot<Unpushed<ReadRequest>>,
    held_event: HeldSlot<Unpushed<SortPhase2Event>>,
    announced: bool,
    sort_stats: bool,
    output_byte_limit: u64,
    /// Fill sizes replacing the resolved ones at `AllAnnounced` (test support).
    #[cfg(test)]
    fills_for_test: Option<(u64, u64)>,
    /// A pool size replacing the resolved one at `AllAnnounced` (test support).
    #[cfg(test)]
    pool_for_test: Option<u64>,
}

impl SpillReadPlanner {
    /// A planner for a merge with `total_memory` (the resolved
    /// `--max-memory`), `phase2_threads` merge-phase threads (outstanding
    /// slices are capped at twice that) and `eligible_clones` `PreadSpillSlices`
    /// clones that can read at once, sharing `supply` with the decompress step
    /// and the merge and reading into `slices` (the pool `PreadSpillSlices`
    /// leases from).
    #[must_use]
    pub fn new(
        total_memory: u64,
        phase2_threads: usize,
        eligible_clones: usize,
        output_byte_limit: u64,
        supply: &SpillSupply,
        slices: Arc<SliceBufferPool>,
    ) -> Self {
        Self {
            slots: Vec::new(),
            by_file_id: HashMap::new(),
            demand: Arc::clone(&supply.demand),
            policy: ReadStreamsPolicy::fixed(1),
            ledger: Arc::clone(&supply.ledger),
            slices,
            total_memory,
            max_inflight_slices: u32::try_from(2 * phase2_threads.max(1)).unwrap_or(u32::MAX),
            cold_inflight_bytes: COLD_INFLIGHT_BYTES,
            eligible_clones: eligible_clones.max(1),
            budget: None,
            hot_prev: Vec::new(),
            cursor: 0,
            next_ordinal: 0,
            outbox: VecDeque::new(),
            held_read: HeldSlot::new(),
            held_event: HeldSlot::new(),
            announced: false,
            sort_stats: false,
            output_byte_limit,
            #[cfg(test)]
            fills_for_test: None,
            #[cfg(test)]
            pool_for_test: None,
        }
    }

    /// Adopt the input's read-stream policy (read only, never observed).
    #[must_use]
    pub fn with_read_streams(mut self, p: Arc<ReadStreamsPolicy>) -> Self {
        self.policy = p;
        self
    }

    /// Log the budget line at `AllAnnounced` (`--sort-stats`).
    #[must_use]
    pub fn with_sort_stats(mut self, on: bool) -> Self {
        self.sort_stats = on;
        self
    }

    /// The cap on outstanding read slices: twice the phase-2 thread count (a
    /// read occupies a pool worker).
    #[must_use]
    pub fn max_inflight_slices(&self) -> u32 {
        self.max_inflight_slices
    }

    /// The resolved budget (test support).
    #[cfg(test)]
    pub(crate) fn budget_for_test(&self) -> Option<ReadAheadBudget> {
        self.budget
    }

    /// Resolve the budget at `AllAnnounced` with these fill sizes (test
    /// support: fills small enough that frames straddle them).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_fills_for_test(mut self, hot_fill: u64, cold_fill: u64) -> Self {
        self.fills_for_test = Some((hot_fill, cold_fill));
        self
    }

    /// Resolve the budget at `AllAnnounced` with this pool size (test
    /// support: make the pool bind).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_pool_for_test(mut self, pool: u64) -> Self {
        self.pool_for_test = Some(pool);
        self
    }

    /// Cap outstanding slices and cold bytes at these values (test support:
    /// force the caps to bind).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_inflight_caps_for_test(mut self, slices: u32, cold_bytes: u64) -> Self {
        self.max_inflight_slices = slices.max(1);
        self.cold_inflight_bytes = cold_bytes;
        self
    }

    /// Run one planning pass and return the requests it queued (test support).
    #[cfg(test)]
    pub(crate) fn plan_pass_for_test(&mut self) -> Vec<ReadRequest> {
        self.plan_pass();
        self.outbox.drain(..).collect()
    }

    /// Handle one phase event (test support).
    #[cfg(test)]
    pub(crate) fn on_event_for_test(&mut self, e: SortPhase1Event) -> SortPhase2Event {
        self.on_event(e)
    }

    /// Register a spill slot at the cold FIFO cap (it holds a hot-sized FIFO
    /// only while in the hot set); an empty one is finalized here without a
    /// read.
    fn register(&mut self, slot: &Arc<SortMergeSlot>) {
        if self.by_file_id.contains_key(&slot.file_id) {
            return;
        }
        slot.set_fifo_cap(FIFO_CAP_COLD);
        let empty = slot.is_empty();
        if empty {
            slot.bp_commit_read(0, true);
            slot.bp_drain_and_finalize();
            self.demand.notify_delivered(slot.file_id);
        }
        self.by_file_id.insert(slot.file_id, self.slots.len());
        self.slots.push(PlannedSlot {
            slot: Arc::clone(slot),
            source: Arc::clone(slot.source()) as Arc<dyn PositionalSource>,
            next_offset: slot.body_start(),
            next_seq: 0,
            eof_issued: empty,
        });
    }

    fn on_event(&mut self, e: SortPhase1Event) -> SortPhase2Event {
        match e {
            SortPhase1Event::SpillReady { slot, path, records_ingested_so_far } => {
                self.register(&slot);
                SortPhase2Event::SpillReady { slot, path, records_ingested_so_far }
            }
            SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far } => {
                SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far }
            }
            SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, total_records } => {
                let budget = ReadAheadBudget::resolve(self.total_memory, slot_count as usize);
                #[cfg(test)]
                let budget = match self.fills_for_test {
                    Some((hot, cold)) => budget.override_for_test(hot, cold),
                    None => budget,
                };
                #[cfg(test)]
                let budget = match self.pool_for_test {
                    Some(pool) => budget.with_pool_for_test(pool),
                    None => budget,
                };
                if self.sort_stats {
                    log::info!("{}", budget.budget_line());
                }
                self.budget = Some(budget);
                self.announced = true;
                SortPhase2Event::AllAnnounced { slot_count, memory_chunk_count, total_records }
            }
        }
    }

    /// The hot set: the merge's demand ([`MergeDemand::hot_ids`]), each slot
    /// only while registered.
    fn hot_set(&self) -> Vec<usize> {
        self.demand
            .hot_ids()
            .as_slice()
            .iter()
            .filter_map(|id| self.by_file_id.get(id))
            .copied()
            .collect()
    }

    /// What slot `s` holds against its allowance: requested bytes plus its
    /// stash charge, less its carried partial frame (which completes only
    /// with the next read, so charging it could block the read that releases
    /// it; the budget counts one carried frame per slot apart).
    fn held(s: &PlannedSlot) -> u64 {
        s.slot.issued_bytes() + s.slot.stash_bytes().saturating_sub(s.slot.carry_bytes())
    }

    /// What slot `s` is charged above its cold allowance (the pool's share).
    fn excess(s: &PlannedSlot, budget: &ReadAheadBudget) -> u64 {
        Self::held(s).saturating_sub(budget.cold_allowance)
    }

    /// What slot `s`'s decompressed FIFO is charged to the pool: its blocks
    /// above the cold cap — the hot cap it was granted, or what a demoted slot
    /// still holds.
    fn fifo_excess(s: &PlannedSlot) -> u64 {
        fifo_charge(u64::from(s.slot.fifo_cap()).max(s.slot.fifo_len_relaxed() as u64))
    }

    /// The pool in use: every slot's charge above its cold terms (read-ahead
    /// and FIFO), plus the idle slice buffers waiting for reuse. One pass over
    /// the slots' mirrors, as the cold rotation already makes.
    fn pool_used(&self, budget: &ReadAheadBudget) -> u64 {
        self.slots.iter().map(|s| Self::excess(s, budget) + Self::fifo_excess(s)).sum::<u64>()
            + self.slices.idle_bytes()
    }

    /// Give hot slot `i` the hot FIFO cap if it lacks it and the pool has room
    /// for the blocks above the cold cap (charging them to `pool_used`); else
    /// it keeps the cold cap, which still lets the merge's awaited slot fill.
    fn grant_hot_fifo(&mut self, i: usize, budget: &ReadAheadBudget, pool_used: &mut u64) {
        let s = &self.slots[i];
        if s.slot.fifo_cap() == FIFO_CAP_HOT {
            return;
        }
        let before = Self::fifo_excess(s);
        let after = fifo_charge(u64::from(FIFO_CAP_HOT).max(s.slot.fifo_len_relaxed() as u64));
        if pool_used.saturating_sub(before) + after <= budget.pool {
            s.slot.set_fifo_cap(FIFO_CAP_HOT);
            *pool_used = pool_used.saturating_sub(before) + after;
        }
    }

    /// One pass: demote the slots that left the hot set, top up the hot set,
    /// then the cold rotation. Returns the number of fills queued.
    fn plan_pass(&mut self) -> u64 {
        let Some(budget) = self.budget else { return 0 };
        let hot = self.hot_set();
        for &i in &self.hot_prev {
            if !hot.contains(&i) {
                self.slots[i].slot.set_fifo_cap(FIFO_CAP_COLD);
            }
        }
        self.hot_prev.clone_from(&hot);
        if let Some(&awaited) = hot.first()
            && self.demand.awaited() == Some(self.slots[awaited].slot.file_id)
        {
            let s = &self.slots[awaited];
            if !s.eof_issued && s.slot.issued_bytes() == 0 && s.slot.stash_len_relaxed() == 0 {
                self.ledger.note_awaited_starved();
            }
        }
        let mut pool_used = self.pool_used(&budget);
        let mut fills = 0;
        for &i in &hot {
            self.grant_hot_fifo(i, &budget, &mut pool_used);
            fills += self.top_up(i, SpillClass::Hot, &budget, &mut pool_used);
        }
        let k = self.slots.len();
        for _ in 0..k {
            if self.ledger.inflight_slices() >= self.max_inflight_slices {
                break;
            }
            let i = self.cursor;
            self.cursor = (self.cursor + 1) % k;
            if !hot.contains(&i) {
                fills += self.top_up(i, SpillClass::Cold, &budget, &mut pool_used);
            }
        }
        fills
    }

    /// The next fill for slot `i` of `class`, or `None` when it may read no
    /// more now. Every slot may read within its cold terms (cold allowance,
    /// cold fills); a hot slot reads past them, in hot fills up to the hot
    /// allowance, only while the pool (`pool_used` of `budget.pool`) has room
    /// for the bytes it would then hold above its cold allowance. A cold slot
    /// already holding more than its cold allowance (a former hot slot) reads
    /// nothing until the merge drains it.
    fn next_fill(
        &self,
        i: usize,
        class: SpillClass,
        budget: &ReadAheadBudget,
        pool_used: u64,
    ) -> Option<u64> {
        let s = &self.slots[i];
        if s.eof_issued || self.ledger.inflight_slices() >= self.max_inflight_slices {
            return None;
        }
        let held = Self::held(s);
        if class == SpillClass::Hot {
            let (hot_allowance, hot_fill) = budget.for_class(SpillClass::Hot);
            let excess = held.saturating_sub(budget.cold_allowance);
            let excess_after = (held + hot_fill).saturating_sub(budget.cold_allowance);
            // `saturating_sub`: the mirrors move under the planner (claims
            // release charges), so this slot's excess can exceed the pass's
            // earlier total.
            if held + hot_fill <= hot_allowance
                && pool_used.saturating_sub(excess) + excess_after <= budget.pool
            {
                return Some(hot_fill);
            }
        } else {
            // The cap never blocks the only outstanding cold read; hot slots
            // skip it, so the awaited slot can always read its cold terms.
            let cold = self.ledger.inflight_cold_bytes();
            if cold != 0 && cold + budget.cold_fill > self.cold_inflight_bytes {
                return None;
            }
        }
        (held + budget.cold_fill <= budget.cold_allowance).then_some(budget.cold_fill)
    }

    /// Keep slot `i` topped up (see [`Self::next_fill`]), charging what it
    /// reads above its cold allowance to `pool_used`. Returns the fills queued.
    fn top_up(
        &mut self,
        i: usize,
        class: SpillClass,
        budget: &ReadAheadBudget,
        pool_used: &mut u64,
    ) -> u64 {
        if class == SpillClass::Cold {
            self.slots[i].slot.set_fifo_cap(FIFO_CAP_COLD);
        }
        let mut fills = 0;
        while let Some(fill) = self.next_fill(i, class, budget, *pool_used) {
            let excess_before = Self::excess(&self.slots[i], budget);
            let s = &mut self.slots[i];
            let len = s.slot.len();
            let want = fill.min(len - s.next_offset);
            let want_usize = usize::try_from(want).expect("a fill fits usize");
            let slices = self.policy.slices_for(want_usize, self.eligible_clones);
            for (offset, slice_len) in slice_ranges(s.next_offset, want_usize, slices) {
                self.outbox.push_back(ReadRequest {
                    ordinal: self.next_ordinal,
                    stream: s.slot.file_id,
                    seq: s.next_seq,
                    source: Arc::clone(&s.source),
                    offset,
                    len: slice_len,
                    last: offset + u64::from(slice_len) == len,
                    target: ReadTarget::Spill { slot: Arc::clone(&s.slot), class },
                });
                self.next_ordinal += 1;
                s.next_seq += 1;
            }
            s.slot.bp_note_issued(want);
            self.ledger.add_inflight(class, want, u32::try_from(slices).unwrap_or(u32::MAX));
            s.next_offset += want;
            if s.next_offset == len {
                s.eof_issued = true;
            }
            *pool_used =
                pool_used.saturating_sub(excess_before) + Self::excess(&self.slots[i], budget);
            fills += 1;
        }
        fills
    }

    /// Whether every registered slot has had its last byte requested.
    fn all_issued(&self) -> bool {
        self.slots.iter().all(|s| s.eof_issued)
    }
}

impl Step for SpillReadPlanner {
    type Input = SortPhase1Event;
    type Outputs = (ReadRequest, SortPhase2Event);

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillReadPlanner",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![
                QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit },
                QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit },
            ],
            branch_ordering: vec![BranchOrdering::None, BranchOrdering::None],
        }
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] =
            &[CounterSpec::new("fills", "fills"), CounterSpec::new("bytes_issued", "bytes")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let view = ctx.outputs.view();
        if let Some(u) = self.held_read.take()
            && let Err(again) = view.a.retry(u)
        {
            self.held_read.put(again);
        }
        if let Some(u) = self.held_event.take()
            && let Err(again) = view.b.retry(u)
        {
            self.held_event.put(again);
        }
        if self.held_read.is_held() || self.held_event.is_held() {
            return Ok(StepOutcome::Contention);
        }
        let mut progressed = false;
        if let Some(e) = ctx.input.pop() {
            let forwarded = self.on_event(e);
            if let Err(u) = view.b.push(forwarded) {
                self.held_event.put(u);
            }
            progressed = true;
        }
        let queued_before = self.outbox.len();
        let fills = self.plan_pass();
        if fills > 0 {
            let fills_bytes: u64 =
                self.outbox.iter().skip(queued_before).map(|r| u64::from(r.len)).sum();
            ctx.counters.add(BYTES_ISSUED, fills_bytes);
            ctx.counters.add(FILLS, fills);
            progressed = true;
        }
        while let Some(req) = self.outbox.pop_front() {
            if let Err(u) = view.a.push(req) {
                self.held_read.put(u);
                break;
            }
            progressed = true;
        }
        if progressed {
            return Ok(StepOutcome::Progress);
        }
        // Retire once every byte is requested. From then on no pass runs, so
        // FIFO caps freeze where they are: a slot hot at that moment keeps the
        // hot cap (still charged to the pool, so the memory bound holds), and
        // one that becomes hot later keeps the cold cap. Only the merge's tail
        // — what is already read ahead — runs under frozen caps; staying alive
        // to keep promoting would cost an O(k) pass per poll for that tail.
        if ctx.input.is_drained()
            && self.announced
            && self.all_issued()
            && self.outbox.is_empty()
            && !self.held_read.is_held()
            && !self.held_event.is_held()
        {
            return Ok(StepOutcome::Finished);
        }
        if ctx.input.is_empty() {
            ctx.input.note_empty_poll();
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests;

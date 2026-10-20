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
//!
//! Each positional fill is advised to the kernel as it is planned
//! ([`PositionalSource::will_need`]): the kernel is asked to start reading the
//! fill's whole range then, so the pool worker that later reads its slices
//! finds some or all of it in the page cache or already in flight, and waits
//! on the device for less. On Linux the advice is issued in 128 KiB calls,
//! because one call starts at most the device's read-ahead window
//! (`fgumi_bam_io::os_hints::advise_willneed_range`). The advice is never
//! binding: whatever the kernel skips or evicts, the read fetches itself.
//!
//! # Consolidation
//!
//! Before `AllAnnounced`, `SpillWrite` may start consolidation jobs
//! (`ConsolidationStart`). The planner opens each job's inputs as slots under
//! fresh slot ids (from [`fgumi_sort::CONSOLIDATION_SLOT_IDS_START`]: a merged
//! run inherits its first input's `file_id` and is read again, so supply maps
//! key by slot id), plans their reads under the job's own budget
//! (`ReadAheadBudget::for_consolidation`), which exists before the merge
//! budget does, and forwards a [`ConsolidationJob`] on its third output to
//! `ConsolidateRuns`. The consolidation merge's demand (awaited, predicted)
//! joins the hot set. On `ConsolidationRetire` the job's slots are marked
//! retired and dropped here, the demands forget them, and the retire epoch
//! moves so the decompress step drops them too; their files close within one
//! pass. At `AllAnnounced` — with no job in flight, or the sort fails — the
//! shared `all_announced` flag is set and the phased cap switches to phase 2
//! before the first merge read is issued, and outstanding slices are capped at
//! twice the phase-2 thread count instead of the phase-1 count.

use std::collections::{HashMap, VecDeque};
use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use fgumi_bam_io::pread::{PositionalSource, ReadStreamsPolicy, SliceBufferPool, slice_ranges};
use fgumi_pipeline_core::{
    CapPhase, HeldSlot, PhasedCap, Unpushed,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};
use fgumi_sort::{MergeDemand, SortMergeSlot};

use super::protocol::{ConsolidationJob, SortPhase1Event, SortPhase2Event};
use super::read_ahead_budget::{COLD_INFLIGHT_BYTES, ReadAheadBudget, fifo_charge};
use super::supply_ledger::{SpillSupply, SupplyLedger};
use crate::pread::{ReadRequest, ReadTarget, SpillClass};

/// Counter slot: fills issued.
const FILLS: usize = 0;
/// Counter slot: bytes issued.
const BYTES_ISSUED: usize = 1;

/// The budget a planned slot reads under.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum BudgetRef {
    /// The merge's, resolved at `AllAnnounced`.
    Merge,
    /// Consolidation job `0`'s.
    Job(u32),
}

/// One registered slot and the planner's cursor into it.
struct PlannedSlot {
    slot: Arc<SortMergeSlot>,
    source: Arc<dyn PositionalSource>,
    next_offset: u64,
    next_seq: u32,
    eof_issued: bool,
    budget: BudgetRef,
    /// The slot's `unissued_total` when the planner last resumed a stream:
    /// what its short reads withheld up to then.
    unissued_seen: u64,
}

/// What one input event forwards.
enum Forward {
    /// A phase event for the merge (branch 1).
    Phase2(SortPhase2Event),
    /// A consolidation job for `ConsolidateRuns` (branch 2).
    Job(ConsolidationJob),
    /// Nothing (a retire).
    None,
}

/// `Serial` planner of the merge's spill reads (see the module docs).
pub struct SpillReadPlanner {
    slots: Vec<PlannedSlot>,
    by_slot_id: HashMap<u32, usize>,
    demand: Arc<MergeDemand>,
    /// The consolidation merge's demand (`None`: no consolidation tail).
    cons_demand: Option<Arc<MergeDemand>>,
    /// Set at `AllAnnounced`; the decompress step reads it to stop notifying
    /// the consolidation demand.
    all_announced: Arc<AtomicBool>,
    /// The supply steps' phased cap, switched to phase 2 at `AllAnnounced`.
    phased: Option<Arc<PhasedCap>>,
    /// Each consolidation job in flight: its budget and its slot ids.
    jobs: HashMap<u32, (ReadAheadBudget, Vec<u32>)>,
    /// The next consolidation slot id to mint.
    next_slot_id: u32,
    /// Bumped on every retire (shared with the decompress step).
    retire_epoch: Arc<AtomicU64>,
    /// Outstanding-slice caps before and after `AllAnnounced`: twice the
    /// phase-1 and twice the phase-2 thread count.
    phase1_inflight_slices: u32,
    phase2_inflight_slices: u32,
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
    held_job: HeldSlot<Unpushed<ConsolidationJob>>,
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
        let inflight = u32::try_from(2 * phase2_threads.max(1)).unwrap_or(u32::MAX);
        Self {
            slots: Vec::new(),
            by_slot_id: HashMap::new(),
            demand: Arc::clone(&supply.demand),
            cons_demand: None,
            all_announced: Arc::new(AtomicBool::new(false)),
            phased: None,
            jobs: HashMap::new(),
            next_slot_id: fgumi_sort::CONSOLIDATION_SLOT_IDS_START,
            retire_epoch: Arc::clone(&supply.retire_epoch),
            phase1_inflight_slices: inflight,
            phase2_inflight_slices: inflight,
            policy: ReadStreamsPolicy::fixed(1),
            ledger: Arc::clone(&supply.ledger),
            slices,
            total_memory,
            max_inflight_slices: inflight,
            cold_inflight_bytes: COLD_INFLIGHT_BYTES,
            eligible_clones: eligible_clones.max(1),
            budget: None,
            hot_prev: Vec::new(),
            cursor: 0,
            next_ordinal: 0,
            outbox: VecDeque::new(),
            held_read: HeldSlot::new(),
            held_event: HeldSlot::new(),
            held_job: HeldSlot::new(),
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

    /// Plan consolidation jobs: `cons_demand` is the consolidation merge's
    /// demand (its awaited and predicted slots join the hot set),
    /// `all_announced` the flag set at `AllAnnounced` (shared with the
    /// decompress step) and `phased` the supply steps' cap, switched to phase 2
    /// there.
    #[must_use]
    pub fn with_consolidation(
        mut self,
        cons_demand: Arc<MergeDemand>,
        all_announced: Arc<AtomicBool>,
        phased: Option<Arc<PhasedCap>>,
    ) -> Self {
        self.cons_demand = Some(cons_demand);
        self.all_announced = all_announced;
        self.phased = phased;
        self
    }

    /// Cap outstanding slices at twice `phase1_threads` until `AllAnnounced`
    /// (the consolidation reads are phase-1 work).
    #[must_use]
    pub fn with_phase1_threads(mut self, phase1_threads: usize) -> Self {
        self.phase1_inflight_slices = u32::try_from(2 * phase1_threads.max(1)).unwrap_or(u32::MAX);
        if !self.announced {
            self.max_inflight_slices = self.phase1_inflight_slices;
        }
        self
    }

    /// Log the budget line at `AllAnnounced` (`--sort-stats`).
    #[must_use]
    pub fn with_sort_stats(mut self, on: bool) -> Self {
        self.sort_stats = on;
        self
    }

    /// The cap on outstanding read slices in the merge: twice the phase-2
    /// thread count (a read occupies a pool worker).
    #[must_use]
    pub fn merge_max_inflight_slices(&self) -> u32 {
        self.phase2_inflight_slices
    }

    /// Every planned slot's id (tests).
    #[cfg(test)]
    pub(crate) fn planned_slot_ids_for_test(&self) -> Vec<u32> {
        self.slots.iter().map(|s| s.slot.slot_id()).collect()
    }

    /// The merge budget, `None` before `AllAnnounced` (tests).
    #[cfg(test)]
    pub(crate) fn merge_budget_for_test(&self) -> Option<ReadAheadBudget> {
        self.budget
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

    /// The hot set's slot ids, in order (test support).
    #[cfg(test)]
    pub(crate) fn hot_set_for_test(&self) -> Vec<u32> {
        self.hot_set().into_iter().map(|i| self.slots[i].slot.slot_id()).collect()
    }

    /// Handle one phase event that forwards a phase-2 event (test support).
    #[cfg(test)]
    pub(crate) fn on_event_for_test(&mut self, e: SortPhase1Event) -> SortPhase2Event {
        match self.on_event(e) {
            Ok(Forward::Phase2(ev)) => ev,
            Ok(Forward::Job(_) | Forward::None) => panic!("not a phase-2 event"),
            Err(e) => panic!("{e}"),
        }
    }

    /// Register a spill slot at the cold FIFO cap (it holds a hot-sized FIFO
    /// only while in the hot set); an empty one is finalized here without a
    /// read.
    fn register(&mut self, slot: &Arc<SortMergeSlot>, budget: BudgetRef, fifo_cold: u32) {
        let id = slot.slot_id();
        if self.by_slot_id.contains_key(&id) {
            return;
        }
        slot.set_fifo_cap(fifo_cold);
        let empty = slot.is_empty();
        if empty {
            slot.bp_commit_read(0, true);
            slot.bp_drain_and_finalize();
            self.notify_delivered(id);
        }
        self.by_slot_id.insert(id, self.slots.len());
        self.slots.push(PlannedSlot {
            slot: Arc::clone(slot),
            source: Arc::clone(slot.source()) as Arc<dyn PositionalSource>,
            next_offset: slot.body_start(),
            next_seq: 0,
            eof_issued: empty,
            budget,
            unissued_seen: 0,
        });
    }

    /// Wake whichever merge awaits `slot_id`.
    fn notify_delivered(&self, slot_id: u32) {
        self.demand.notify_delivered(slot_id);
        if let Some(d) = &self.cons_demand {
            d.notify_delivered(slot_id);
        }
    }

    /// The budget slot `s` reads under, if it exists yet.
    fn budget_of(&self, s: &PlannedSlot) -> Option<ReadAheadBudget> {
        match s.budget {
            BudgetRef::Merge => self.budget,
            BudgetRef::Job(id) => self.jobs.get(&id).map(|(b, _)| *b),
        }
    }

    /// Open a job's inputs as slots under fresh slot ids, register them under
    /// the job's budget and build the job for `ConsolidateRuns`.
    fn start_job(
        &mut self,
        job_id: u32,
        inputs: Vec<super::protocol::ConsolidationInput>,
        key_kind: fgumi_sort::SpillKeyKind,
    ) -> io::Result<ConsolidationJob> {
        if self.announced {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("consolidation job {job_id} started after AllAnnounced"),
            ));
        }
        let budget = ReadAheadBudget::for_consolidation(self.total_memory, inputs.len());
        self.ledger.note_consolidation_job(&budget);
        let mut slots = Vec::with_capacity(inputs.len());
        let mut input_paths = Vec::with_capacity(inputs.len());
        for input in inputs {
            let slot_id = self.next_slot_id;
            self.next_slot_id = slot_id.checked_add(1).ok_or_else(|| {
                io::Error::other("consolidation slot ids exhausted (more than 2^31 inputs)")
            })?;
            // Opening reads the file's magic: an external wait on the
            // filesystem (#1068).
            let slot = fgumi_pipeline_core::external_wait(|| {
                fgumi_sort::open_spill_slot_as(&input.path, input.file_id, slot_id)
            })
            .map_err(|e| {
                io::Error::other(format!(
                    "SpillReadPlanner: failed to open consolidation input {}: {e:#}",
                    input.path.display()
                ))
            })?;
            self.register(&slot, BudgetRef::Job(job_id), budget.fifo_cold);
            slots.push(slot);
            input_paths.push(input.path);
        }
        let ids = slots.iter().map(|s| s.slot_id()).collect();
        self.jobs.insert(job_id, (budget, ids));
        Ok(ConsolidationJob { job_id, slots, key_kind, input_paths })
    }

    /// Retire a completed job's slots: mark them, drop them here, make the
    /// demands forget them and move the retire epoch for the decompress step.
    fn retire_job(&mut self, job_id: u32) -> io::Result<()> {
        let Some((_, ids)) = self.jobs.remove(&job_id) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("ConsolidationRetire for unknown consolidation job {job_id}"),
            ));
        };
        for s in &self.slots {
            if s.budget == BudgetRef::Job(job_id) {
                s.slot.retire();
            }
        }
        for &id in &ids {
            self.demand.forget(id);
            if let Some(d) = &self.cons_demand {
                d.forget(id);
            }
        }
        self.slots.retain(|s| s.budget != BudgetRef::Job(job_id));
        self.by_slot_id =
            self.slots.iter().enumerate().map(|(i, s)| (s.slot.slot_id(), i)).collect();
        self.hot_prev.clear();
        self.cursor = 0;
        self.retire_epoch.fetch_add(1, Ordering::Release);
        Ok(())
    }

    fn on_event(&mut self, e: SortPhase1Event) -> io::Result<Forward> {
        Ok(match e {
            SortPhase1Event::SpillReady { slot, path, records_ingested_so_far } => {
                self.register(&slot, BudgetRef::Merge, super::read_ahead_budget::FIFO_CAP_COLD);
                Forward::Phase2(SortPhase2Event::SpillReady { slot, path, records_ingested_so_far })
            }
            SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far } => {
                Forward::Phase2(SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far })
            }
            SortPhase1Event::ConsolidationStart { job_id, inputs, key_kind } => {
                Forward::Job(self.start_job(job_id, inputs, key_kind)?)
            }
            SortPhase1Event::ConsolidationRetire { job_id } => {
                self.retire_job(job_id)?;
                Forward::None
            }
            SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, total_records } => {
                if let Some(&id) = self.jobs.keys().min() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("AllAnnounced arrived with consolidation job {id} in flight"),
                    ));
                }
                // Before the first merge read is planned: every phase-2 read,
                // decompress and compress runs under the phase-2 cap.
                self.all_announced.store(true, Ordering::Release);
                if let Some(p) = &self.phased {
                    p.switch_to(CapPhase::Two);
                }
                self.max_inflight_slices = self.phase2_inflight_slices;
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
                Forward::Phase2(SortPhase2Event::AllAnnounced {
                    slot_count,
                    memory_chunk_count,
                    total_records,
                })
            }
        })
    }

    /// The hot set: the merge's demand ([`MergeDemand::hot_ids`]: awaited,
    /// predicted and frontier slots, deduplicated), each only while registered
    /// and not yet EOF (at most [`super::read_ahead_budget::HOT_SLOTS`]).
    fn hot_set(&self) -> Vec<usize> {
        let merge = self.demand.hot_ids();
        let cons = self.cons_demand.as_ref().map(|d| [d.awaited(), d.predicted()]);
        let mut hot = Vec::new();
        let ids = merge.as_slice().iter().copied().chain(cons.into_iter().flatten().flatten());
        for id in ids {
            if let Some(&i) = self.by_slot_id.get(&id)
                && !hot.contains(&i)
                && !self.slots[i].slot.queue_eof()
            {
                hot.push(i);
            }
        }
        hot
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
    fn pool_used(&self, group: BudgetRef, budget: &ReadAheadBudget) -> u64 {
        self.slots
            .iter()
            .filter(|s| s.budget == group)
            .map(|s| Self::excess(s, budget) + Self::fifo_excess(s))
            .sum::<u64>()
            + self.slices.idle_bytes()
    }

    /// Give hot slot `i` the hot FIFO cap if it lacks it and the pool has room
    /// for the blocks above the cold cap (charging them to `pool_used`); else
    /// it keeps the cold cap, which still lets the merge's awaited slot fill.
    fn grant_hot_fifo(&mut self, i: usize, budget: &ReadAheadBudget, pool_used: &mut u64) {
        let s = &self.slots[i];
        if s.slot.fifo_cap() == budget.fifo_hot {
            return;
        }
        let before = Self::fifo_excess(s);
        let after = fifo_charge(u64::from(budget.fifo_hot).max(s.slot.fifo_len_relaxed() as u64));
        if pool_used.saturating_sub(before) + after <= budget.pool {
            s.slot.set_fifo_cap(budget.fifo_hot);
            *pool_used = pool_used.saturating_sub(before) + after;
        }
    }

    /// One pass: demote the slots that left the hot set, top up the hot set,
    /// then the cold rotation, each slot under its own budget (the merge's, or
    /// its consolidation job's; one not yet resolved plans nothing). Returns
    /// the number of fills queued.
    fn plan_pass(&mut self) -> u64 {
        // A stream's slot leaves the table once its EOF was ingested: it will
        // never be read again.
        if self.slots.iter().any(|s| s.slot.is_sequential() && s.slot.eof_seen()) {
            for s in &mut self.slots {
                if s.slot.is_sequential() && s.slot.eof_seen() {
                    s.eof_issued = true;
                }
            }
            self.slots.retain(|s| !(s.slot.is_sequential() && s.slot.eof_seen()));
            self.by_slot_id =
                self.slots.iter().enumerate().map(|(i, s)| (s.slot.slot_id(), i)).collect();
            self.hot_prev.clear();
            self.cursor = 0;
        }
        if self.slots.is_empty() {
            return 0;
        }
        let hot = self.hot_set();
        for &i in &self.hot_prev {
            if !hot.contains(&i)
                && let Some(b) = self.budget_of(&self.slots[i])
            {
                self.slots[i].slot.set_fifo_cap(b.fifo_cold);
            }
        }
        self.hot_prev.clone_from(&hot);
        if let Some(&awaited) = hot.first()
            && self.demand.awaited() == Some(self.slots[awaited].slot.slot_id())
        {
            let s = &self.slots[awaited];
            if !s.eof_issued && s.slot.issued_bytes() == 0 && s.slot.stash_len_relaxed() == 0 {
                self.ledger.note_awaited_starved();
            }
        }
        // Each budget's pool in use, computed once per pass on first touch.
        let mut pool_used: HashMap<BudgetRef, u64> = HashMap::new();
        let mut fills = 0;
        for &i in &hot {
            let group = self.slots[i].budget;
            let Some(budget) = self.budget_of(&self.slots[i]) else { continue };
            let mut u = match pool_used.get(&group) {
                Some(&u) => u,
                None => self.pool_used(group, &budget),
            };
            self.grant_hot_fifo(i, &budget, &mut u);
            fills += self.top_up(i, Self::class_of(group, true), &budget, &mut u);
            pool_used.insert(group, u);
        }
        let k = self.slots.len();
        for _ in 0..k {
            if self.ledger.inflight_slices() >= self.max_inflight_slices {
                break;
            }
            let i = self.cursor;
            self.cursor = (self.cursor + 1) % k;
            if hot.contains(&i) {
                continue;
            }
            let group = self.slots[i].budget;
            let Some(budget) = self.budget_of(&self.slots[i]) else { continue };
            let mut u = match pool_used.get(&group) {
                Some(&u) => u,
                None => self.pool_used(group, &budget),
            };
            fills += self.top_up(i, Self::class_of(group, false), &budget, &mut u);
            pool_used.insert(group, u);
        }
        fills
    }

    /// The read class of a slot in `group`, hot or not.
    fn class_of(group: BudgetRef, hot: bool) -> SpillClass {
        match (group, hot) {
            (BudgetRef::Merge, true) => SpillClass::Hot,
            (BudgetRef::Merge, false) => SpillClass::Cold,
            (BudgetRef::Job(_), true) => SpillClass::ConsolidationHot,
            (BudgetRef::Job(_), false) => SpillClass::ConsolidationCold,
        }
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
        if class.is_hot() {
            let (hot_allowance, hot_fill) = budget.for_class(class);
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
        if !class.is_hot() {
            self.slots[i].slot.set_fifo_cap(budget.fifo_cold);
        }
        let mut fills = 0;
        let sequential = self.slots[i].slot.is_sequential();
        if sequential {
            // A stream is read in order, one read in flight: the next read is
            // issued only once the previous one is ingested, and none after a
            // zero read (the stream's EOF).
            let s = &mut self.slots[i];
            match s.slot.sequential_read_gate(s.next_seq) {
                fgumi_sort::SequentialReadGate::Eof => {
                    s.eof_issued = true;
                    return 0;
                }
                fgumi_sort::SequentialReadGate::InFlight => return 0,
                fgumi_sort::SequentialReadGate::Due => {}
            }
            // A read returns what the stream has, so it may end short of its
            // fill: resume at the bytes it delivered, which is where the
            // stream is. Every slice is ingested (above), so the total is
            // exact.
            let unissued = s.slot.unissued_total();
            s.next_offset -= unissued - s.unissued_seen;
            s.unissued_seen = unissued;
            // And only once the stream has bytes (or has ended): the read then
            // takes what is there without waiting, so no pool worker waits on
            // a producer, which may itself be waiting for another stream's
            // read to be issued to a worker (a coupled producer, such as
            // `tee`). A stream with nothing yet is planned again next pass:
            // its pump wakes no one when bytes arrive, so this relies on the
            // runtime polling an idle planner (its workers' parks all time
            // out), as it already does for the merge's demand.
            if !s.source.ready_at(s.next_offset) {
                return 0;
            }
        }
        while let Some(fill) = self.next_fill(i, class, budget, *pool_used) {
            let excess_before = Self::excess(&self.slots[i], budget);
            let s = &mut self.slots[i];
            let len = s.slot.len();
            let want = fill.min(len - s.next_offset);
            let want_usize = usize::try_from(want).expect("a fill fits usize");
            let slices = if sequential {
                1
            } else {
                // Ask the kernel to start reading the whole fill now, so the
                // pool worker that later reads each slice finds it in the page
                // cache or in flight instead of starting the device read
                // itself. Bounded to the fill just planned; a stream has
                // nothing to advise. Issued here, at plan time, for every
                // fill (hot ones included: they are the awaited ones), so the
                // device read gets the lead time the queued request waits.
                //
                // The advised pages are page cache: they are charged to no
                // budget and the kernel may evict them before the slice is
                // read. That is acceptable because the hint is advisory only;
                // a miss degrades to the blocking read there was before.
                //
                // Cost, on the serial planner: measured warm on macOS (one
                // `F_RDADVISE` per fill), a three-input 5.3 GB merge advised
                // 1,264 fills in 44-46 ms (~35 us a fill) and a small
                // consolidating sort 5.9-6.4 k fills in 146-165 ms (~25 us),
                // on a planner that stayed under 0.5% busy; wall time with the
                // hint was no worse than without it in either run, two runs
                // each. That is not material, so there is no size threshold.
                // Linux advises a 4 MiB fill in 32 calls; their cost on Linux
                // was not measured.
                s.source.will_need(s.next_offset, want);
                self.policy.slices_for(want_usize, self.eligible_clones)
            };
            for (offset, slice_len) in slice_ranges(s.next_offset, want_usize, slices) {
                self.outbox.push_back(ReadRequest {
                    ordinal: self.next_ordinal,
                    stream: s.slot.slot_id(),
                    seq: s.next_seq,
                    source: Arc::clone(&s.source),
                    offset,
                    len: slice_len,
                    last: !sequential && offset + u64::from(slice_len) == len,
                    sequential,
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
            if sequential {
                break;
            }
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
    type Outputs = (ReadRequest, SortPhase2Event, ConsolidationJob);

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillReadPlanner",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![
                QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit },
                QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit },
                QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit },
            ],
            branch_ordering: vec![BranchOrdering::None, BranchOrdering::None, BranchOrdering::None],
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
        if let Some(u) = self.held_job.take()
            && let Err(again) = view.c.retry(u)
        {
            self.held_job.put(again);
        }
        if self.held_read.is_held() || self.held_event.is_held() || self.held_job.is_held() {
            return Ok(StepOutcome::Contention);
        }
        let mut progressed = false;
        if let Some(e) = ctx.input.pop() {
            match self.on_event(e)? {
                Forward::Phase2(ev) => {
                    if let Err(u) = view.b.push(ev) {
                        self.held_event.put(u);
                    }
                }
                Forward::Job(job) => {
                    if let Err(u) = view.c.push(job) {
                        self.held_job.put(u);
                    }
                }
                Forward::None => {}
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
            && !self.held_job.is_held()
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

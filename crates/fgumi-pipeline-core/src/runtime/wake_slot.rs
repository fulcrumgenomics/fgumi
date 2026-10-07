//! Pipeline thread identity for wake routing and hold timing.
//!
//! Every pipeline OS thread has a *wake slot*: pool worker `w` is slot `w`,
//! detached driver `d` is slot `n_threads + d` (the numbering the
//! `WorkerStateBoard` already uses). `Pipeline::run` enters the slot on the
//! spawned thread; a queue's reject path reads it to record *which* thread is
//! holding an item. Threads outside a pipeline loop have no slot.

use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

// Under `--cfg loom` the words, flags and fences of the wake protocols below
// (`HolderSet`, `DirectParked`, `ThreadSlots`, `delivery_fence`) are loom's,
// and so is the thread handle a slot holds, so `tests/loom_wake.rs`
// model-checks the real code: record → fence → re-check, arm → fence →
// re-poll, register → fence → first pass, and the producer's fence before it
// reads any of them. `HoldClock` stays on std atomics: it is statistics, not
// part of any wake protocol.
#[cfg(loom)]
use loom::sync::atomic::{AtomicBool as SlotFlag, AtomicU64 as HolderWord, fence as holder_fence};
#[cfg(not(loom))]
use std::sync::atomic::{AtomicU64 as HolderWord, fence as holder_fence};

/// The OS-thread handle a directed wake unparks (`loom`'s under `--cfg loom`).
#[cfg(not(loom))]
pub type ParkThread = std::thread::Thread;
/// The OS-thread handle a directed wake unparks (`loom`'s under `--cfg loom`).
#[cfg(loom)]
pub type ParkThread = loom::thread::Thread;

/// The calling thread's [`ParkThread`] handle.
#[must_use]
pub fn current_thread() -> ParkThread {
    #[cfg(not(loom))]
    return std::thread::current();
    #[cfg(loom)]
    return loom::thread::current();
}

/// The producer half of the unpark-path fence rule: issued after a publish
/// (a push, a pop that made room, a limit raise) and before reading who to
/// wake (a [`ThreadSlots`] slot, a [`HolderSet`]). Pairs with the fence in
/// [`ThreadSlots::register`], [`HolderSet::latch_and_recheck`] and
/// [`DirectParked::arm`]: of each pair, whichever runs second sees the other's
/// write.
#[inline]
pub fn delivery_fence() {
    holder_fence(Ordering::SeqCst);
}

use crate::padded::Padded;
use crate::runtime::stats::PipelineStats;
use crate::topology::StepIdx;

thread_local! {
    static SLOT: Cell<Option<usize>> = const { Cell::new(None) };
}

/// Sets this OS thread's wake slot for its lifetime; clears it on drop.
pub(crate) struct SlotGuard(());

impl SlotGuard {
    /// Enter wake slot `slot` on the calling thread.
    pub(crate) fn enter(slot: usize) -> Self {
        SLOT.with(|s| s.set(Some(slot)));
        Self(())
    }
}

impl Drop for SlotGuard {
    fn drop(&mut self) {
        SLOT.with(|s| s.set(None));
        clear_thread_flags();
    }
}

/// Clear every per-thread wake hint ([`take_rejected`], [`take_flushed`],
/// [`take_cap_refused`], [`take_resume_step`], the cap record masks), so a hint set by one pipeline
/// loop on a thread cannot leak into the next loop that thread runs (the
/// caller's thread runs worker 0 of a one-worker run, or the fused path).
/// Called when a thread leaves its slot and after the fused loop.
pub(crate) fn clear_thread_flags() {
    POPPED.with(|p| p.set(false));
    EC_ARMED.with(|a| a.set(false));
    REJECTED.with(|r| r.set(false));
    FLUSHED.with(|f| f.set(false));
    CAP_REFUSED.with(|r| r.set(false));
    RESUME_STEP.with(|r| r.set(None));
    CAP_RECORDED.with(|r| r.set(0));
    CAP_RECORDED_ANY.with(|r| r.set(false));
    CAP_OWED.with(|r| r.set(0));
}

/// This OS thread's wake slot, if it is a pipeline thread.
#[inline]
pub(crate) fn current_slot() -> Option<usize> {
    SLOT.with(Cell::get)
}

thread_local! {
    static POPPED: Cell<bool> = const { Cell::new(false) };
}

/// The current dispatch took an input item (`BranchInputHandle::pop` returned
/// one). Set on every successful pop; taken by `dispatch_one_step` after the
/// dispatch.
#[inline]
pub(crate) fn note_popped() {
    POPPED.with(|p| p.set(true));
}

/// Whether this thread popped an input item since the last call; clears it.
/// A dispatch that took input reports `Progress` (the `StepOutcome::Progress`
/// contract): the dispatcher upgrades an idle outcome after a pop, so the
/// reverse wake that releases an upstream holder for the freed slot is never
/// lost to a step that forgot.
#[inline]
pub(crate) fn take_popped() -> bool {
    POPPED.with(|p| p.replace(false))
}

thread_local! {
    static REJECTED: Cell<bool> = const { Cell::new(false) };
    static FLUSHED: Cell<bool> = const { Cell::new(false) };
}

/// A tracked edge (one with a holder set) refused this thread's push for good,
/// so the thread is holding an item and must idle where a direct unpark
/// reaches it. Sets the thread-local flag [`take_rejected`] reads. Called only
/// when the edge has a holder set: a Legacy byte-bounded edge with only a hold
/// clock (`--pipeline-stats`) must leave the idle branch on the event-count.
#[inline]
pub(crate) fn note_held() {
    REJECTED.with(|r| r.set(true));
}

/// Whether an edge with a holder set refused a push by this thread for good
/// ([`note_held`]) since the last call. The worker loop reads it when it goes
/// idle: a thread that may be holding an item must idle where a direct
/// `unpark` reaches it.
#[inline]
pub(crate) fn take_rejected() -> bool {
    REJECTED.with(|r| r.replace(false))
}

/// Non-consuming peek at the flag [`take_rejected`] takes, for tests (same
/// shape as [`cap_refused_pending`]).
#[cfg(test)]
pub(crate) fn rejected_pending() -> bool {
    REJECTED.with(Cell::get)
}

/// A held item was just re-pushed successfully (`BranchOutputHandle::retry`
/// returned `Ok`) by this thread. `dispatch_one_step` takes it after the
/// dispatch, so in a Directed plan the forward wake of that push runs even if
/// the step then reports `NoProgress`, `Contention` or `Capped` (a Legacy plan
/// ignores it).
#[inline]
pub(crate) fn note_flushed() {
    FLUSHED.with(|f| f.set(true));
}

/// Whether this thread flushed a held item since the last call; clears it.
/// Same shape as [`take_rejected`] and [`take_cap_refused`]: one thread-local
/// load and store.
#[inline]
pub(crate) fn take_flushed() -> bool {
    FLUSHED.with(|f| f.replace(false))
}

/// Non-consuming peek, for tests (same shape as [`cap_refused_pending`]).
#[cfg(test)]
pub(crate) fn flushed_pending() -> bool {
    FLUSHED.with(Cell::get)
}

thread_local! {
    static EC_ARMED: Cell<bool> = const { Cell::new(false) };
}

/// This thread is (`true`) or no longer is (`false`) armed on the pool
/// event-count: between its `prepare_wait` and the `cancel_wait`/`wait` that
/// balances it. Set by the worker loop around the arm→re-poll→wait window.
#[inline]
pub(crate) fn set_ec_armed(armed: bool) {
    EC_ARMED.with(|a| a.set(armed));
}

/// Whether this thread is armed on the pool event-count (see
/// [`set_ec_armed`]): a `request_worker` made from a step in its re-poll must
/// not count the requester itself as a parked worker to wake.
#[inline]
pub(crate) fn ec_armed() -> bool {
    EC_ARMED.with(Cell::get)
}

thread_local! {
    static CAP_REFUSED: Cell<bool> = const { Cell::new(false) };
    static RESUME_STEP: Cell<Option<StepIdx>> = const { Cell::new(None) };
}

/// A phase cap refused this thread and recorded its wake slot: the cap's next
/// release will `unpark` it (`PhaseCap::release`). Set by the cap's refusal
/// path, read by the worker loop when it goes idle.
#[inline]
pub(crate) fn note_cap_refused() {
    CAP_REFUSED.with(|r| r.set(true));
}

/// Whether a phase cap recorded this thread since the last call. A worker
/// whose last pass was refused this way must idle where the cap release's
/// direct `unpark` reaches it — on its timer, not the event-count.
#[inline]
pub(crate) fn take_cap_refused() -> bool {
    CAP_REFUSED.with(|r| r.replace(false))
}

/// Whether a phase cap recorded this thread since the flag was last taken,
/// without clearing it.
#[inline]
pub(crate) fn cap_refused_pending() -> bool {
    CAP_REFUSED.with(Cell::get)
}

thread_local! {
    static CAP_RECORDED: Cell<u64> = const { Cell::new(0) };
    static CAP_RECORDED_ANY: Cell<bool> = const { Cell::new(false) };
    static CAP_OWED: Cell<u64> = const { Cell::new(0) };
}

/// A record of this thread now stands on some cap (a refusal or a full-cap
/// skip, or a skip whose fresh record a release claimed). Sets the "recorded
/// anything" flag, which a thread on the shared anonymous bit also needs (it
/// must park like any recorded thread), and — when the record is the
/// thread's own bit on the wake plan's cap `bit` — bit `bit` of the
/// per-thread mask the next pass start takes (`WakePlan::start_pass`). So a
/// non-zero mask implies the flag.
#[inline]
pub(crate) fn note_cap_recorded(bit: Option<u32>) {
    CAP_RECORDED_ANY.with(|r| r.set(true));
    if let Some(i) = bit {
        CAP_RECORDED.with(|r| r.set(r.get() | (1 << i)));
    }
}

/// The caps this thread recorded on since the mask was last taken; clears it,
/// and the "recorded anything" flag with it.
#[inline]
pub(crate) fn take_cap_recorded() -> u64 {
    CAP_RECORDED_ANY.with(|r| r.set(false));
    CAP_RECORDED.with(|r| r.replace(0))
}

/// Whether this thread recorded on any cap since the mask was last taken,
/// with or without a bit in the mask ([`note_cap_recorded`] sets the flag
/// for both).
#[inline]
pub(crate) fn cap_recorded_pending() -> bool {
    CAP_RECORDED_ANY.with(Cell::get)
}

/// Set the caps whose wake this thread owes for the pass now starting: it
/// was recorded there, and a release claimed the record (woke this thread).
#[inline]
pub(crate) fn set_cap_owed(mask: u64) {
    CAP_OWED.with(|r| r.set(mask));
}

/// This thread was admitted on the wake plan's cap `i`: it spends a permit
/// there, so a wake it owes that cap is paid (its own release wakes the next
/// recorded thread). One thread-local read; a write only when owed.
#[inline]
pub(crate) fn note_cap_admitted(i: u32) {
    CAP_OWED.with(|r| {
        let owed = r.get();
        if owed & (1 << i) != 0 {
            r.set(owed & !(1 << i));
        }
    });
}

/// The caps whose wake this thread still owes at the end of a pass; clears it.
#[inline]
pub(crate) fn take_cap_owed() -> u64 {
    CAP_OWED.with(|r| r.replace(0))
}

/// Remember the capped step that refused this thread, so the thread's next
/// pass polls it before anything else: a release wakes one refused worker per
/// freed permit, and that worker must spend the permit, not round-robin to
/// other work and leave it idle.
#[inline]
pub(crate) fn set_resume_step(step: StepIdx) {
    RESUME_STEP.with(|r| r.set(Some(step)));
}

/// Take the step [`set_resume_step`] recorded, if any.
#[inline]
pub(crate) fn take_resume_step() -> Option<StepIdx> {
    RESUME_STEP.with(Cell::take)
}

/// The threads waiting on one resource, as one bit per wake slot (plus one
/// shared, anonymous bit for threads without a slot). Two users:
///
/// - an edge's holders: set by a refused producer ([`Self::latch_and_recheck`],
///   or [`Self::record_current`] under the caller's lock); cleared by
///   [`Self::take`] (the consumer's `WakePlan::on_progress`, the rebalancer's
///   `WakePlan::wake_holders`), which unparks what it clears, so a producer
///   whose retry succeeds leaves at most one spurious wake behind;
/// - a phase cap's refused workers (`PhaseCap`): set by a refusal's or a
///   full-cap skip's latch; cleared by `take_up_to` (a release, or a
///   forwarded wake, claims one per freed permit and wakes it), by
///   [`Self::clear_slot`] (an admitted thread withdraws the record its own
///   call set, with no wake), by [`Self::clear_slot_reporting`] (the owner's
///   pass start) and by [`Self::clear_all`] (a new run binds the cap).
///
/// The anonymous bit is only ever cleared by a claim, which broadcasts.
pub struct HolderSet {
    /// One cache line per word: a refusal or a claim writes one word, and the
    /// words of different edges and caps must not share a line.
    words: Box<[Padded<HolderWord>]>,
    n_slots: usize,
}

impl HolderSet {
    /// A set for `n_slots` pipeline threads plus the shared anonymous slot.
    #[must_use]
    pub fn new(n_slots: usize) -> Self {
        Self {
            words: (0..(n_slots + 1).div_ceil(64)).map(|_| Padded(HolderWord::new(0))).collect(),
            n_slots,
        }
    }

    /// The slot recorded for a thread without a wake slot.
    #[cfg(test)]
    #[must_use]
    pub fn anonymous_slot(&self) -> usize {
        self.n_slots
    }

    /// Record `slot`'s bit (`Relaxed`). A slot outside the set's range, or a
    /// thread with no slot, is recorded in the anonymous bit, which the waker
    /// broadcasts (`WakePlan::deliver_anonymous`).
    #[inline]
    fn record_slot(&self, slot: Option<usize>) -> bool {
        let bit = self.bit_of(slot);
        let mask = 1u64 << (bit % 64);
        self.words[bit / 64].0.fetch_or(mask, Ordering::Relaxed) & mask == 0
    }

    /// Record the calling thread with no fence, and without marking the thread.
    /// For a refusal decided under the caller's own lock
    /// (`ItemQueue::note_rejected_holder`, the reorder stash cap): that lock,
    /// not a fence, orders this store before the consumer's take, because the
    /// in-order pop takes the same lock first. The caller calls [`note_held`]
    /// itself.
    #[inline]
    pub(crate) fn record_current(&self) {
        self.record_slot(current_slot());
    }

    /// A push into this set's edge was just refused. Record the calling thread,
    /// fence, and re-check: `true` means room appeared between the refusal and
    /// the record, and the caller retries instead of holding. The `SeqCst`
    /// fence pairs with the consumer's fence before `take`
    /// (`WakePlan::on_progress`) or the rebalancer's (`WakePlan::wake_holders`):
    /// of the two, whichever runs second sees the other's write. The bit is not
    /// cleared here: on an edge only `take` clears it (the consumer's
    /// `on_progress`, or the rebalancer's `wake_holders`), so a retry that
    /// succeeds leaves at most one spurious unpark.
    #[inline]
    pub(crate) fn latch_and_recheck(&self, has_room: impl FnOnce() -> bool) -> bool {
        self.latch_slot(current_slot(), has_room)
    }

    /// [`Self::latch_and_recheck`] with an explicit slot, for the loom models,
    /// whose threads never enter a pipeline slot (as
    /// `PhaseCap::try_acquire_as`), and for the unit test that pins the
    /// record → predicate order.
    #[cfg(any(test, loom))]
    #[doc(hidden)]
    pub fn latch_as(&self, slot: usize, has_room: impl FnOnce() -> bool) -> bool {
        self.latch_slot(Some(slot), has_room)
    }

    /// The whole reject protocol, in its one required order: record, then
    /// fence, then evaluate the predicate. Evaluating the predicate before the
    /// record reopens the lost wake (a take between the two finds no bit, and
    /// the predicate missed the room), so every transport, and every phase cap
    /// (`PhaseCap`, which records a slot it is handed), goes through here.
    #[inline]
    pub(crate) fn latch_slot(&self, slot: Option<usize>, has_room: impl FnOnce() -> bool) -> bool {
        self.latch_slot_reporting(slot, has_room).1
    }

    /// [`Self::latch_slot`], also reporting whether this record set the bit
    /// (`false`: it was already set, by an earlier record still standing).
    #[inline]
    pub(crate) fn latch_slot_reporting(
        &self,
        slot: Option<usize>,
        has_room: impl FnOnce() -> bool,
    ) -> (bool, bool) {
        #[cfg(test)]
        crate::queues::test_hooks::run_before_latch();
        let newly_set = self.record_slot(slot);
        holder_fence(Ordering::SeqCst);
        #[cfg(test)]
        crate::queues::test_hooks::run_after_record();
        (newly_set, has_room())
    }

    /// The bit `slot` is recorded in: its own, or the anonymous one for a
    /// slot outside the set's range or a thread with none.
    #[inline]
    fn bit_of(&self, slot: Option<usize>) -> usize {
        slot.filter(|&s| s < self.n_slots).unwrap_or(self.n_slots)
    }

    /// Clear `slot`'s own bit (the recorded thread took what it waited for on
    /// its own). The anonymous bit is shared, so it is never cleared here.
    #[inline]
    pub(crate) fn clear_slot(&self, slot: Option<usize>) {
        let bit = self.bit_of(slot);
        if bit != self.n_slots {
            self.words[bit / 64].0.fetch_and(!(1u64 << (bit % 64)), Ordering::Relaxed);
        }
    }

    /// Whether `slot`'s bit (or the anonymous one it maps to) is set (tests).
    #[cfg(test)]
    pub(crate) fn is_set(&self, slot: Option<usize>) -> bool {
        let bit = self.bit_of(slot);
        self.words[bit / 64].0.load(Ordering::Relaxed) & (1u64 << (bit % 64)) != 0
    }

    /// Whether `slot` has a bit of its own (not the shared anonymous one), so
    /// its owner can tell whether a release claimed it.
    #[inline]
    pub(crate) fn owns_bit(&self, slot: Option<usize>) -> bool {
        self.bit_of(slot) != self.n_slots
    }

    /// [`Self::clear_slot`] after a `Relaxed` check that the bit is set (so a
    /// thread with no record pays a load, not a read-modify-write), reporting
    /// whether this call cleared it.
    #[inline]
    pub(crate) fn clear_slot_reporting(&self, slot: Option<usize>) -> bool {
        let bit = self.bit_of(slot);
        if bit == self.n_slots {
            return false;
        }
        let mask = 1u64 << (bit % 64);
        self.words[bit / 64].0.load(Ordering::Relaxed) & mask != 0
            && self.words[bit / 64].0.fetch_and(!mask, Ordering::Relaxed) & mask != 0
    }

    /// Clear every bit (a new run binds the set afresh).
    pub(crate) fn clear_all(&self) {
        for w in &self.words {
            w.0.store(0, Ordering::Relaxed);
        }
    }

    /// Whether no bit is recorded: one `Relaxed` load per word.
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.words.iter().all(|w| w.0.load(Ordering::Relaxed) == 0)
    }

    /// Whether `slot` is the anonymous bit (a broadcast, not a thread).
    #[inline]
    pub(crate) fn is_anonymous(&self, slot: usize) -> bool {
        slot == self.n_slots
    }

    /// Claim recorded holders bit by bit in ascending order, calling `f(slot)`
    /// for each, until `n` are claimed; returns how many. A bit another
    /// claimer cleared first is skipped. One `Relaxed` load per word when
    /// nothing is recorded. The caller fences (`SeqCst`) first.
    pub(crate) fn take_up_to(&self, n: usize, f: &mut dyn FnMut(usize)) -> usize {
        claim_bits(&self.words, None, None, n, &mut |slot| {
            f(slot);
            true
        })
    }

    /// Take and clear every recorded holder, calling `f(slot)` for each in
    /// ascending order; returns how many. One `Relaxed` load per word; the
    /// `swap` runs only on a non-zero word. The caller fences (`SeqCst`) first.
    #[inline]
    pub fn take(&self, f: &mut dyn FnMut(usize)) -> usize {
        let mut n = 0;
        for (i, w) in self.words.iter().enumerate() {
            if w.0.load(Ordering::Relaxed) == 0 {
                continue;
            }
            let mut bits = w.0.swap(0, Ordering::Relaxed);
            while bits != 0 {
                let b = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                f(i * 64 + b);
                n += 1;
            }
        }
        n
    }
}

/// One registered OS-thread handle per directed-wake target (a plan's drivers,
/// or its pool workers), empty until that thread registers itself.
///
/// The registration protocol: the thread registers its own handle, fences
/// (`SeqCst`) and only then starts its first pass, while a producer publishes,
/// fences ([`delivery_fence`]) and then reads the slot. So a publish either
/// finds the handle and unparks it, or is seen by the first pass. A handle
/// registered by another thread after `spawn` has no such ordering against the
/// thread's first pass.
pub struct ThreadSlots {
    slots: Box<[ThreadSlot]>,
}

/// A slot: `std`'s `OnceLock` (its `get` is an `Acquire` load of the state its
/// `set` releases), or, under loom, a loom flag with the same orderings in
/// front of a `OnceLock` loom does not see.
#[derive(Default)]
struct ThreadSlot {
    #[cfg(loom)]
    set: SlotFlag,
    thread: std::sync::OnceLock<ParkThread>,
}

impl ThreadSlot {
    fn get(&self) -> Option<&ParkThread> {
        #[cfg(loom)]
        if !self.set.load(Ordering::Acquire) {
            return None;
        }
        self.thread.get()
    }

    fn set(&self, t: ParkThread) {
        let _ = self.thread.set(t);
        #[cfg(loom)]
        self.set.store(true, Ordering::Release);
    }
}

impl ThreadSlots {
    /// `n` empty slots.
    #[must_use]
    pub fn new(n: usize) -> Self {
        Self { slots: (0..n).map(|_| ThreadSlot::default()).collect() }
    }

    /// How many slots.
    #[must_use]
    pub(crate) fn len(&self) -> usize {
        self.slots.len()
    }

    /// Register `t` in slot `i` (the first registration wins; out of range is
    /// ignored), then fence. Called by the registering thread itself, before
    /// its first pass: the fence orders the registration before that pass.
    pub fn register(&self, i: usize, t: ParkThread) {
        if let Some(slot) = self.slots.get(i) {
            slot.set(t);
            holder_fence(Ordering::SeqCst);
        }
    }

    /// Unpark the thread registered in slot `i`; `false` (nothing done) when
    /// the slot is empty or out of range. The caller has fenced
    /// ([`delivery_fence`]) after its publish.
    #[inline]
    pub fn unpark(&self, i: usize) -> bool {
        match self.slots.get(i).and_then(ThreadSlot::get) {
            Some(t) => {
                t.unpark();
                true
            }
            None => false,
        }
    }

    /// Unpark every registered thread.
    pub fn unpark_all(&self) {
        for slot in &self.slots {
            if let Some(t) = slot.get() {
                t.unpark();
            }
        }
    }

    /// Whether slot `i` holds a registered thread.
    #[cfg(test)]
    #[must_use]
    pub fn is_registered(&self, i: usize) -> bool {
        self.slots.get(i).and_then(ThreadSlot::get).is_some()
    }
}

/// The pool workers parked on their own timer that a `Pool` wake may unpark
/// when no event-count waiter exists: one bit per worker.
///
/// The direct-park protocol: the worker arms its bit, fences (`SeqCst`) and
/// re-polls before it parks; the producer publishes, fences (inside
/// `PoolEventCount::notify_one`), finds no event-count waiter and claims an
/// armed bit. So a publish either reaches an armed worker or is seen by its
/// re-poll.
pub struct DirectParked {
    /// One cache line per word, as `HolderSet`'s.
    words: Box<[Padded<HolderWord>]>,
}

impl DirectParked {
    /// No worker armed, for `n_workers` workers.
    #[must_use]
    pub fn new(n_workers: usize) -> Self {
        Self { words: (0..n_workers.div_ceil(64)).map(|_| Padded(HolderWord::new(0))).collect() }
    }

    /// Arm worker `w` (it is about to re-poll, then park on its timer), then
    /// fence: the fence orders the arm before the re-poll. `false` (nothing
    /// armed) when `w` is out of range.
    pub fn arm(&self, w: usize) -> bool {
        let Some(word) = self.words.get(w / 64) else { return false };
        word.0.fetch_or(1u64 << (w % 64), Ordering::Relaxed);
        holder_fence(Ordering::SeqCst);
        true
    }

    /// A claim mask for [`Self::claim_one`] with one bit set per worker in
    /// `workers`, in the armed words' layout (bit `w % 64` of word `w / 64`)
    /// for `n_workers` workers.
    ///
    /// # Panics
    ///
    /// Panics if a worker in `workers` is `>= n_workers`.
    #[must_use]
    pub fn mask_for(workers: &[usize], n_workers: usize) -> Box<[u64]> {
        let mut mask = vec![0u64; n_workers.div_ceil(64)];
        for &w in workers {
            assert!(w < n_workers, "worker {w} out of range ({n_workers} workers)");
            mask[w / 64] |= 1u64 << (w % 64);
        }
        mask.into_boxed_slice()
    }

    /// Worker `w` is running again.
    pub fn disarm(&self, w: usize) {
        if let Some(word) = self.words.get(w / 64) {
            word.0.fetch_and(!(1u64 << (w % 64)), Ordering::Relaxed);
        }
    }

    /// Whether worker `w` is armed (tests).
    #[cfg(test)]
    pub(crate) fn is_armed(&self, w: usize) -> bool {
        self.words
            .get(w / 64)
            .is_some_and(|word| word.0.load(Ordering::Relaxed) & (1 << (w % 64)) != 0)
    }

    /// Whether no worker is armed: one `Relaxed` load per word.
    #[inline]
    #[must_use]
    pub(crate) fn is_empty(&self) -> bool {
        self.words.iter().all(|w| w.0.load(Ordering::Relaxed) == 0)
    }

    /// Claim armed workers in ascending order, clearing each bit, until
    /// `unpark(w)` returns `true` for one (it reached a thread); `false` when
    /// none did. `mask` (one bit per worker, as the armed words; a missing word
    /// is all-zero) restricts the claim to its set bits; `None`: any armed
    /// worker. The caller has fenced after its publish.
    pub fn claim_one(&self, mask: Option<&[u64]>, unpark: &mut dyn FnMut(usize) -> bool) -> bool {
        claim_bits(&self.words, mask, None, 1, unpark) == 1
    }

    /// [`Self::claim_one`] over every armed worker, leaving worker `skip`'s bit
    /// armed and unclaimed (the caller itself, which must not wake itself).
    pub fn claim_one_except(
        &self,
        skip: Option<usize>,
        unpark: &mut dyn FnMut(usize) -> bool,
    ) -> bool {
        claim_bits(&self.words, None, skip, 1, unpark) == 1
    }
}

/// The one claim loop of the wake bit sets: claim set bits in ascending order,
/// clearing each with `fetch_and` (a bit another claimer cleared first is
/// skipped; bit `skip`, if any, is left alone), and call `f(bit)` on each
/// claimed one; stop once `limit` calls returned `true`. Returns how many did.
/// `eligible` (`None`: every bit) restricts the claim to its set bits. One
/// `Relaxed` load per word when nothing is set. The caller fences (`SeqCst`)
/// first.
fn claim_bits(
    words: &[Padded<HolderWord>],
    eligible: Option<&[u64]>,
    skip: Option<usize>,
    limit: usize,
    f: &mut dyn FnMut(usize) -> bool,
) -> usize {
    let mut counted = 0;
    if limit == 0 {
        return 0;
    }
    for (i, w) in words.iter().enumerate() {
        let allowed = eligible.map_or(u64::MAX, |m| m.get(i).copied().unwrap_or(0));
        let mut bits = w.0.load(Ordering::Relaxed) & allowed;
        if let Some(s) = skip.filter(|&s| s / 64 == i) {
            bits &= !(1u64 << (s % 64));
        }
        while bits != 0 {
            let b = bits.trailing_zeros() as usize;
            let mask = 1u64 << b;
            bits &= !mask;
            if w.0.fetch_and(!mask, Ordering::Relaxed) & mask != 0 && f(i * 64 + b) {
                counted += 1;
                if counted == limit {
                    return counted;
                }
            }
        }
    }
    counted
}

/// Times held items on one byte-bounded edge (stats-on only). It keeps the
/// time of the first rejection of the hold in progress; the next successful
/// push into the edge (or a stash insert that accepts the item) ends the hold,
/// and rejections in between are retries.
///
/// The scope of a hold depends on how the producer runs. A `Parallel` step has
/// one clone per pool worker, each holding its own item and retrying it on its
/// own thread, so the clock keeps one stamp per thread slot
/// ([`Self::per_thread`]). Every other producer holds at most one item per
/// output branch, and a `Serial` step's held item may be retried by whichever
/// worker next wins the step's lock, so the clock keeps one stamp for the step
/// ([`Self::per_step`]): a flush by another worker ends the hold. One padded
/// word per stamp, so producers on different threads never share a line here.
pub struct HoldClock {
    step: StepIdx,
    stats: Arc<PipelineStats>,
    origin: Instant,
    /// Nanoseconds since `origin` (+1, so 0 means "not holding"). Per thread
    /// slot plus a last, shared entry for threads that have none; or a single
    /// entry for the whole step.
    since: Box<[Padded<AtomicU64>]>,
}

impl HoldClock {
    /// A clock that books holds to the `Parallel` producer `step` in `stats`,
    /// with one stamp per pipeline thread (`n_slots`) plus one shared stamp for
    /// threads without a slot.
    #[must_use]
    pub fn per_thread(step: StepIdx, stats: Arc<PipelineStats>, n_slots: usize) -> Self {
        Self::with_stamps(step, stats, n_slots + 1)
    }

    /// A clock that books holds to producer `step` in `stats` with one stamp
    /// for the step, whichever thread rejects or pushes: for every producer
    /// that is not `Parallel`.
    #[must_use]
    pub fn per_step(step: StepIdx, stats: Arc<PipelineStats>) -> Self {
        Self::with_stamps(step, stats, 1)
    }

    fn with_stamps(step: StepIdx, stats: Arc<PipelineStats>, n: usize) -> Self {
        Self {
            step,
            stats,
            origin: Instant::now(),
            since: (0..n).map(|_| Padded(AtomicU64::new(0))).collect(),
        }
    }

    /// How many stamps the clock keeps: `n_slots + 1` per thread, 1 per step.
    #[cfg(test)]
    pub(crate) fn stamps(&self) -> usize {
        self.since.len()
    }

    /// The stamp of the calling thread (per-thread clock) or of the step.
    #[inline]
    fn slot(&self) -> &AtomicU64 {
        let last = self.since.len() - 1;
        if last == 0 {
            return &self.since[0].0;
        }
        let s = current_slot().filter(|&s| s < last).unwrap_or(last);
        &self.since[s].0
    }

    #[inline]
    fn now_ns(&self) -> u64 {
        u64::try_from(self.origin.elapsed().as_nanos()).unwrap_or(u64::MAX).saturating_add(1)
    }

    /// A push was rejected: the first rejection of a hold stamps it, later ones
    /// count a retry.
    #[inline]
    pub(crate) fn on_reject(&self) {
        let c = self.slot();
        if c.load(Ordering::Relaxed) == 0 {
            c.store(self.now_ns(), Ordering::Relaxed);
        } else {
            self.stats.record_held_retry(self.step);
        }
    }

    /// A push succeeded (or a stash accepted the item): if a hold was open, it
    /// ends now.
    #[inline]
    pub(crate) fn on_push(&self) {
        let c = self.slot();
        let t = c.load(Ordering::Relaxed);
        if t != 0 {
            c.store(0, Ordering::Relaxed);
            self.stats.record_hold(self.step, self.now_ns().saturating_sub(t));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slot_guard_sets_and_clears_the_thread_slot() {
        assert_eq!(current_slot(), None);
        let g = SlotGuard::enter(3);
        assert_eq!(current_slot(), Some(3));
        let other = std::thread::spawn(current_slot).join().unwrap();
        assert_eq!(other, None, "another thread never sees this thread's slot");
        drop(g);
        assert_eq!(current_slot(), None);
    }

    #[test]
    fn holder_set_take_visits_each_bit_once() {
        let set = std::sync::Arc::new(HolderSet::new(130));
        assert_eq!(set.anonymous_slot(), 130);
        let mut handles: Vec<_> = [0usize, 63, 64, 129]
            .into_iter()
            .map(|slot| {
                let set = std::sync::Arc::clone(&set);
                std::thread::spawn(move || {
                    let _g = SlotGuard::enter(slot);
                    set.record_current();
                })
            })
            .collect();
        let set_anon = std::sync::Arc::clone(&set);
        handles.push(std::thread::spawn(move || set_anon.record_current()));
        for h in handles {
            h.join().unwrap();
        }
        let mut seen = Vec::new();
        assert_eq!(set.take(&mut |s| seen.push(s)), 5);
        assert_eq!(seen, vec![0, 63, 64, 129, 130]);
        let mut again = Vec::new();
        assert_eq!(set.take(&mut |s| again.push(s)), 0);
        assert!(again.is_empty());
    }

    /// The reject protocol's order, pinned deterministically: the predicate
    /// runs after the record (it takes, and must find its own bit), a retry
    /// leaves its bit for the next take, and an out-of-range slot is recorded
    /// as the anonymous bit. `record_current` marks nothing: only `note_held`
    /// does.
    #[test]
    fn latch_as_records_before_it_evaluates_the_predicate() {
        let set = HolderSet::new(2);
        let mut in_predicate = Vec::new();
        let retry = set.latch_as(1, || {
            set.take(&mut |s| in_predicate.push(s));
            false
        });
        assert!(!retry, "no room: hold");
        assert_eq!(in_predicate, vec![1], "recorded before the predicate ran");
        assert!(set.latch_as(7, || true), "room: retry");
        let mut seen = Vec::new();
        set.take(&mut |s| seen.push(s));
        assert_eq!(seen, vec![set.anonymous_slot()]);
        set.record_current();
        assert!(!rejected_pending(), "a record alone does not mark the thread");
    }

    /// A claimed worker with no thread to unpark is skipped: the claim goes on
    /// to the next armed worker, so one `Pool` wake is not lost on it.
    #[test]
    fn a_claim_moves_past_an_armed_worker_it_cannot_unpark() {
        let parked = DirectParked::new(2);
        assert!(parked.arm(0) && parked.arm(1));
        let mut tried = Vec::new();
        assert!(parked.claim_one(None, &mut |w| {
            tried.push(w);
            w == 1
        }));
        assert_eq!(tried, vec![0, 1], "worker 0 had no thread; worker 1 was woken");
        assert!(!parked.claim_one(None, &mut |_| true), "both bits were claimed");
    }

    /// `mask_for` sets exactly the named workers' bits, across words.
    #[test]
    fn mask_for_sets_one_bit_per_worker_across_words() {
        let mask = DirectParked::mask_for(&[1, 64, 65, 129], 130);
        assert_eq!(&*mask, &[1u64 << 1, 0b11, 1u64 << 1]);
        assert_eq!(&*DirectParked::mask_for(&[], 130), &[0, 0, 0]);
    }

    /// A masked claim crosses words and claims only masked bits: with workers
    /// 0, 63, 64 and 129 armed and the mask {63, 129}, the claims reach 63 then
    /// 129, and 0 and 64 stay armed. A mask shorter than the armed words reads
    /// its missing words as zero.
    #[test]
    fn a_masked_claim_spans_words_and_claims_only_masked_workers() {
        let parked = DirectParked::new(130);
        for w in [0, 63, 64, 129] {
            assert!(parked.arm(w));
        }
        let mask = DirectParked::mask_for(&[63, 129], 130);
        let mut woke = Vec::new();
        assert!(parked.claim_one(Some(&mask), &mut |w| {
            woke.push(w);
            true
        }));
        assert!(parked.claim_one(Some(&mask), &mut |w| {
            woke.push(w);
            true
        }));
        assert!(!parked.claim_one(Some(&mask), &mut |_| true), "no masked worker is left");
        assert_eq!(woke, vec![63, 129]);
        assert!(parked.is_armed(0) && parked.is_armed(64), "unmasked workers stay armed");

        // One word of mask over three armed words: only word 0 is eligible.
        assert!(parked.arm(129));
        let short = [1u64];
        let mut woke = Vec::new();
        assert!(parked.claim_one(Some(&short), &mut |w| {
            woke.push(w);
            true
        }));
        assert_eq!(woke, vec![0]);
        assert!(!parked.claim_one(Some(&short), &mut |_| true), "word 2 is outside the mask");
        assert!(parked.is_armed(64) && parked.is_armed(129));
    }

    /// Leaving a slot clears every per-thread hint, so none leaks into the
    /// thread's next pipeline loop.
    #[test]
    fn leaving_a_slot_clears_every_thread_hint() {
        let g = SlotGuard::enter(1);
        note_held();
        note_flushed();
        note_cap_refused();
        note_popped();
        note_cap_recorded(Some(4));
        set_cap_owed(1 << 4);
        set_resume_step(StepIdx(3));
        drop(g);
        assert!(!rejected_pending() && !flushed_pending() && !cap_refused_pending());
        assert!(!take_popped());
        assert!(!cap_recorded_pending(), "the record mask and flag");
        assert_eq!(take_cap_owed(), 0, "the owed mask");
        assert_eq!(take_resume_step(), None);
    }

    #[test]
    fn flushed_is_taken_once() {
        let _ = take_flushed();
        note_flushed();
        assert!(flushed_pending());
        assert!(take_flushed());
        assert!(!take_flushed());
        assert!(!flushed_pending());
    }
}

//! Merge-wide demand facts shared by the k-way merge consumer (`SortMerge` /
//! `MergeDriver`) and the spill supply (`SortSpillDecompress`): which file the
//! merge is waiting on and the consumer's thread handle, so the worker that
//! delivers the awaited block can unpark it. Only a delivery to the awaited
//! file wakes the merge: waking on every delivery costs a park per unrelated
//! block (v0.7.0 paid ~1.9 parks per stall that way). Created once per sort by
//! the chain builder's `add_sort`.
//!
//! Lost-wakeup argument: [`MergeDemand::await_slot`] stores `awaited` and THEN
//! re-checks the slot under its `decompressed` mutex; the producer pushes under
//! that mutex and THEN reads `awaited`. Mutex release/acquire orders the two
//! critical sections, so either the consumer's re-check sees the block or the
//! producer sees the awaited id (and `Thread::unpark` leaves a token if the
//! park has not started). `await_slot` is the only consumer-side sequence; the
//! loom model `loom_merge_wake_never_lost` in `tests/loom_merge_slots.rs` calls
//! it, so swapping its store and its re-check fails the model.

use std::sync::atomic::AtomicBool as StdAtomicBool;
use std::sync::atomic::AtomicU64 as StdAtomicU64;
use std::sync::atomic::Ordering::Relaxed;

use crossbeam_utils::CachePadded;
#[cfg(loom)]
use loom::sync::Mutex;
#[cfg(loom)]
use loom::sync::atomic::{AtomicU64, Ordering};
#[cfg(loom)]
use loom::thread::{self as thread_mod, Thread};
#[cfg(not(loom))]
use std::sync::Mutex;
#[cfg(not(loom))]
use std::sync::atomic::{AtomicU64, Ordering};
#[cfg(not(loom))]
use std::thread::{self as thread_mod, Thread};

use crate::merge_slots::SortMergeSlot;

/// What the awaited slot was doing when the merge stalled on it
/// ([`SortMergeSlot::awaited_state`]). Not `merge_stalls`' finer-grained
/// stall census state, which classifies the owned engine's awaited file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AwaitedSlotState {
    /// No block is being decompressed for the slot (`in_flight == 0`).
    Starved,
    /// Blocks are being decompressed for the slot (`in_flight > 0`).
    Decompressing,
}

impl AwaitedSlotState {
    /// The number of states (the stats' bucket count).
    const COUNT: usize = 2;
}

/// Outcome of the merge's pool-worker request, mirrored from pipeline-core's
/// `PoolRequest`: fgumi-sort does not depend on pipeline-core, and these
/// counters must exist whenever `--sort-stats` does. The same outcomes are
/// also counted per step by pipeline-core (the `req_*` wake-table columns),
/// deliberately: that table needs pipeline stats on and has no
/// `Unavailable` column, and neither can relate a request to the stall it was
/// made in, which the sleeper share here does.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MergePoolRequest {
    /// A parked worker was woken.
    Woken,
    /// Every worker was awake.
    AllAwake,
    /// An unacknowledged request was outstanding.
    Pending,
    /// No pool to ask.
    Unavailable,
}

impl MergePoolRequest {
    /// The number of outcomes (the stats' bucket count).
    const COUNT: usize = 4;
}

/// The merge's demand counters (relaxed; read once at the end). Plain std
/// atomics even under loom: they are not part of the wake protocol. Written by
/// the consumer, except `wakes_delivered`, which a producer bumps only when it
/// actually unparks the consumer.
#[derive(Default)]
pub struct MergeDemandStats {
    stall_episodes: StdAtomicU64,
    stall_ns: StdAtomicU64,
    registrations: StdAtomicU64,
    parking_registrations: StdAtomicU64,
    wakes_delivered: StdAtomicU64,
    /// Indexed by `AwaitedSlotState as usize`.
    awaited: [StdAtomicU64; AwaitedSlotState::COUNT],
    /// Indexed by `MergePoolRequest as usize`.
    pool: [StdAtomicU64; MergePoolRequest::COUNT],
    requests_with_sleeper: StdAtomicU64,
    stall_ns_with_sleeper: StdAtomicU64,
    partial_flushes: StdAtomicU64,
    /// The supply decompresses inline (`--sort::file-granularity`), so stalls
    /// are not classified: see [`Self::mark_inline_decompress`].
    inline_decompress: StdAtomicBool,
}

impl MergeDemandStats {
    /// One stall episode began with the awaited slot in `state`.
    pub fn record_stall(&self, state: AwaitedSlotState) {
        self.stall_episodes.fetch_add(1, Relaxed);
        self.awaited[state as usize].fetch_add(1, Relaxed);
    }

    /// A stall episode ended after `ns` nanoseconds.
    pub fn record_stall_ns(&self, ns: u64) {
        self.stall_ns.fetch_add(ns, Relaxed);
    }

    /// The consumer parked after its last registration (its `await_slot`
    /// returned `false` and the dispatch delivered nothing, so its driver
    /// parks until the awaited slot's delivery or its idle timer).
    pub fn record_park(&self) {
        self.parking_registrations.fetch_add(1, Relaxed);
    }

    /// The merge asked the pool for a worker and got `r`, with
    /// `parked_workers` workers parked on the pool's event-count when it
    /// asked. Returns whether the request found a sleeping worker: one parked
    /// on the event-count, or — `r == Woken` with none there — one woken from
    /// a timer park (a pinned or holding worker, a cap-parked one, or worker 0
    /// at one thread), which the event-count census does not see.
    pub fn record_pool_request(&self, r: MergePoolRequest, parked_workers: usize) -> bool {
        let sleeper = parked_workers > 0 || r == MergePoolRequest::Woken;
        self.pool[r as usize].fetch_add(1, Relaxed);
        if sleeper {
            self.requests_with_sleeper.fetch_add(1, Relaxed);
        }
        sleeper
    }

    /// A stall episode whose pool request found a sleeping worker lasted `ns`.
    pub fn record_stall_with_sleeper_ns(&self, ns: u64) {
        self.stall_ns_with_sleeper.fetch_add(ns, Relaxed);
    }

    /// The merge flushed a partial output batch because it stalled.
    pub fn record_partial_flush(&self) {
        self.partial_flushes.fetch_add(1, Relaxed);
    }

    /// The merge's supply decompresses inline, one worker per file
    /// (`--sort::file-granularity`). That path never raises a slot's
    /// `in_flight`, and the merge classifies a stall from lock-free slot state
    /// only (it must not take the reader lock), so every stall reports
    /// `Starved`. The awaited-slot line says so instead of printing a
    /// misleading 100% starved.
    pub fn mark_inline_decompress(&self) {
        self.inline_decompress.store(true, Relaxed);
    }
}

/// A point-in-time copy of [`MergeDemandStats`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MergeDemandSnapshot {
    /// Stall episodes (one per pull that stalled before producing a record).
    pub stall_episodes: u64,
    /// Total time spent in stall episodes, in nanoseconds.
    pub stall_ns: u64,
    /// Times the consumer declared a slot awaited.
    pub registrations: u64,
    /// Registrations after which the consumer parked: its `await_slot` found
    /// nothing and the dispatch delivered nothing. A registration whose
    /// re-check found a block or EOF, or whose stall flushed a partial batch
    /// (the merge runs again at once), is not one of these.
    /// `wakes_delivered` is measured against this count.
    pub parking_registrations: u64,
    /// Times a producer's delivery unparked the consumer.
    pub wakes_delivered: u64,
    /// Stall episodes whose awaited slot had nothing being decompressed.
    pub awaited_starved: u64,
    /// Stall episodes whose awaited slot had blocks being decompressed.
    pub awaited_decompressing: u64,
    /// Pool requests that woke a parked worker.
    pub pool_woken: u64,
    /// Pool requests that found every worker awake.
    pub pool_all_awake: u64,
    /// Pool requests refused because one was already pending.
    pub pool_pending: u64,
    /// Pool requests with no pool to ask.
    pub pool_unavailable: u64,
    /// Pool requests that found a sleeping worker: parked on the pool's
    /// event-count when the merge asked, or woken from a timer park.
    pub requests_with_sleeper: u64,
    /// Stall time of episodes whose request found a sleeping worker, in ns.
    pub stall_ns_with_sleeper: u64,
    /// Partial output batches flushed on a stall.
    pub partial_flushes: u64,
    /// The supply decompresses inline, so the awaited-slot buckets are not
    /// classified (every stall reports `Starved`).
    pub inline_decompress: bool,
}

/// `num / den` as a whole percentage, rounded to nearest; `0` when `den == 0`.
fn pct_int(num: u64, den: u64) -> u64 {
    if den == 0 { 0 } else { (num * 100 + den / 2) / den }
}

/// `num / den` as a percentage; `0.0` when `den == 0`.
#[allow(clippy::cast_precision_loss)]
fn pct(num: u64, den: u64) -> f64 {
    if den == 0 { 0.0 } else { num as f64 * 100.0 / den as f64 }
}

/// Nanoseconds as seconds.
#[allow(clippy::cast_precision_loss)]
fn secs(ns: u64) -> f64 {
    ns as f64 / 1e9
}

impl MergeDemandSnapshot {
    /// The `--sort-stats` lines, in their fixed order: demand, awaited slot
    /// state, pool requests at stalls, partial flushes. Line 2's percentages
    /// are of `stall_episodes`; line 3's sleeper share is of pool requests.
    #[must_use]
    pub fn log_lines(&self) -> Vec<String> {
        let requests =
            self.pool_woken + self.pool_all_awake + self.pool_pending + self.pool_unavailable;
        let n = self.stall_episodes;
        vec![
            format!(
                "Merge demand: {} stall episodes ({:.1} s exact), wakes delivered {} of {} \
                 parking registrations ({} registrations)",
                n,
                secs(self.stall_ns),
                self.wakes_delivered,
                self.parking_registrations,
                self.registrations
            ),
            format!(
                "Awaited slot at stall: starved {}% / decompressing {}%{}",
                pct_int(self.awaited_starved, n),
                pct_int(self.awaited_decompressing, n),
                if self.inline_decompress {
                    " (inline decompress: stalls are not classified, every one reports starved)"
                } else {
                    ""
                }
            ),
            format!(
                "Pool at stall: requests woken {} / all-awake {} / pending {} / unavailable {}; \
                 a worker was asleep at {:.1}% of requests ({:.1} s)",
                self.pool_woken,
                self.pool_all_awake,
                self.pool_pending,
                self.pool_unavailable,
                pct(self.requests_with_sleeper, requests),
                secs(self.stall_ns_with_sleeper)
            ),
            format!("Merge output: {} partial flushes on stall", self.partial_flushes),
        ]
    }
}

/// Merge-wide demand state; see the module doc.
pub struct MergeDemand {
    /// `file_id + 1`; [`NONE`] = the consumer is not waiting. On its own cache
    /// line: every decompress worker reads it on every delivery, and nothing
    /// else that is written often may share its line.
    awaited: CachePadded<AtomicU64>,
    /// The thread to unpark when the awaited file's block lands.
    consumer: Mutex<Option<Thread>>,
    /// Raised while the merge is parked on an empty awaited slot; the pool
    /// scheduler reads it (the sort's refill hint) to run the spill supply
    /// before other pool work. A plain std atomic even under loom: it steers
    /// scheduling only and is not part of the wake protocol.
    starved: std::sync::Arc<std::sync::atomic::AtomicBool>,
    stats: MergeDemandStats,
}

/// `awaited` value meaning "not waiting".
const NONE: u64 = 0;

impl Default for MergeDemand {
    fn default() -> Self {
        Self::new()
    }
}

impl MergeDemand {
    /// Empty demand: nothing awaited, no consumer registered.
    #[must_use]
    pub fn new() -> Self {
        Self {
            awaited: CachePadded::new(AtomicU64::new(NONE)),
            consumer: Mutex::new(None),
            starved: std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false)),
            stats: MergeDemandStats::default(),
        }
    }

    /// The shared `starved` flag, for the pool scheduler's refill hint.
    #[must_use]
    pub fn starved_signal(&self) -> std::sync::Arc<std::sync::atomic::AtomicBool> {
        std::sync::Arc::clone(&self.starved)
    }

    /// Raise or lower `starved`. Stores only on a change, so a merge that
    /// re-stalls within one episode does not rewrite a line every pool worker
    /// reads on every pass.
    pub fn set_starved(&self, on: bool) {
        if self.starved.load(Relaxed) != on {
            self.starved.store(on, Relaxed);
        }
    }

    /// Record the calling thread as the consumer to unpark. Cheap when already
    /// registered (one `ThreadId` compare under an uncontended lock); replaces
    /// the stored handle when the merge runs on a different thread, so a stale
    /// handle can never swallow a wake. Half of [`Self::await_slot`]; crate
    /// private so no caller can register without the re-check.
    ///
    /// # Panics
    ///
    /// Panics if the consumer mutex is poisoned.
    pub(crate) fn register_consumer(&self) {
        let me = thread_mod::current();
        let mut slot = self.consumer.lock().expect("consumer mutex poisoned");
        if slot.as_ref().map(Thread::id) != Some(me.id()) {
            *slot = Some(me);
        }
    }

    /// Declare that the merge waits on `file_id`. Half of
    /// [`Self::await_slot`], which re-checks the slot afterwards: a bare
    /// `set_awaited` followed by a park can lose a wake, so it is crate
    /// private.
    pub(crate) fn set_awaited(&self, file_id: u32) {
        self.awaited.store(u64::from(file_id) + 1, Ordering::Release);
        self.stats.registrations.fetch_add(1, Relaxed);
    }

    /// The merge is no longer waiting: a later delivery does not unpark it.
    /// Safe at any time — the merge re-registers through [`Self::await_slot`]
    /// (and so re-checks) before it parks on a slot again.
    pub fn clear_awaited(&self) {
        self.awaited.store(NONE, Ordering::Release);
    }

    /// The awaited file, if any.
    ///
    /// # Panics
    ///
    /// Never in practice: `awaited` only ever holds `NONE` or a `u32` id + 1.
    #[must_use]
    pub fn awaited(&self) -> Option<u32> {
        match self.awaited.load(Ordering::Acquire) {
            NONE => None,
            v => Some(u32::try_from(v - 1).expect("awaited encodes a u32 file id")),
        }
    }

    /// The consumer's side of the wake protocol: register this thread, declare
    /// `slot` awaited, then re-check it under its `decompressed` mutex. Returns
    /// `true` — and clears `awaited` — when the slot can already make progress
    /// (the caller must not park); `false` when the caller may park, because any
    /// later delivery to `slot` will see `awaited` and unpark it. Every call
    /// counts a registration; the caller records a park
    /// ([`MergeDemandStats::record_park`]) when it does park.
    ///
    /// # Panics
    ///
    /// Panics if the consumer mutex or the slot's `decompressed` mutex is
    /// poisoned.
    pub fn await_slot(&self, slot: &SortMergeSlot) -> bool {
        self.register_consumer();
        self.set_awaited(slot.file_id);
        if slot.has_block_or_eof() {
            self.clear_awaited();
            return true;
        }
        false
    }

    /// A producer delivered a block, EOF or failure to `file_id` (called AFTER
    /// releasing the slot's `decompressed` mutex). Unparks the consumer iff
    /// `file_id` is the awaited file and a consumer is registered; `true` iff a
    /// thread was unparked. The common case — another file — is one `Acquire`
    /// load on a line nobody else writes; the CAS runs only on a match, so a
    /// delivery to an unawaited file never takes the line exclusive. The load
    /// runs after the producer released `decompressed`, where the lost-wakeup
    /// argument (module doc) needs it.
    ///
    /// # Panics
    ///
    /// Panics if the consumer mutex is poisoned.
    pub fn notify_delivered(&self, file_id: u32) -> bool {
        let id = u64::from(file_id) + 1;
        if self.awaited.load(Ordering::Acquire) != id {
            return false;
        }
        if self.awaited.compare_exchange(id, NONE, Ordering::AcqRel, Ordering::Acquire).is_err() {
            return false;
        }
        let woke = match self.consumer.lock().expect("consumer mutex poisoned").as_ref() {
            Some(t) => {
                t.unpark();
                true
            }
            None => false,
        };
        if woke {
            self.stats.wakes_delivered.fetch_add(1, Relaxed);
        }
        woke
    }

    /// The demand counters.
    #[must_use]
    pub fn stats(&self) -> &MergeDemandStats {
        &self.stats
    }

    /// A point-in-time copy of every counter.
    #[must_use]
    pub fn snapshot(&self) -> MergeDemandSnapshot {
        let s = &self.stats;
        let ld = |a: &StdAtomicU64| a.load(Relaxed);
        MergeDemandSnapshot {
            stall_episodes: ld(&s.stall_episodes),
            stall_ns: ld(&s.stall_ns),
            registrations: ld(&s.registrations),
            parking_registrations: ld(&s.parking_registrations),
            wakes_delivered: ld(&s.wakes_delivered),
            awaited_starved: ld(&s.awaited[AwaitedSlotState::Starved as usize]),
            awaited_decompressing: ld(&s.awaited[AwaitedSlotState::Decompressing as usize]),
            pool_woken: ld(&s.pool[MergePoolRequest::Woken as usize]),
            pool_all_awake: ld(&s.pool[MergePoolRequest::AllAwake as usize]),
            pool_pending: ld(&s.pool[MergePoolRequest::Pending as usize]),
            pool_unavailable: ld(&s.pool[MergePoolRequest::Unavailable as usize]),
            requests_with_sleeper: ld(&s.requests_with_sleeper),
            stall_ns_with_sleeper: ld(&s.stall_ns_with_sleeper),
            partial_flushes: ld(&s.partial_flushes),
            inline_decompress: s.inline_decompress.load(Relaxed),
        }
    }

    /// The registered consumer's thread id (tests only).
    #[cfg(all(test, not(loom)))]
    pub(crate) fn consumer_id_for_test(&self) -> Option<std::thread::ThreadId> {
        self.consumer.lock().unwrap().as_ref().map(Thread::id)
    }
}

// Outside `loom::model` the loom primitives are illegal, so these compile only
// in a normal build (as `merge_slots`' tests do).
#[cfg(all(test, not(loom)))]
mod tests {
    use std::sync::Arc;
    use std::sync::mpsc::channel;
    use std::time::Duration;

    use super::*;

    #[test]
    fn notify_for_a_non_awaited_file_does_not_wake() {
        let d = MergeDemand::new();
        d.register_consumer();
        d.set_awaited(3);
        assert!(!d.notify_delivered(4));
        assert_eq!(d.awaited(), Some(3));
        assert_eq!(d.snapshot().wakes_delivered, 0);
    }

    #[test]
    fn notify_for_the_awaited_file_wakes_exactly_once() {
        let d = MergeDemand::new();
        d.register_consumer();
        d.set_awaited(7);
        assert!(d.notify_delivered(7));
        assert!(!d.notify_delivered(7), "the CAS cleared it; a second delivery does not wake");
        assert_eq!(d.awaited(), None);
        assert_eq!(d.snapshot().wakes_delivered, 1);
        assert_eq!(d.snapshot().registrations, 1);
    }

    /// With no consumer registered there is nobody to unpark: the awaited id is
    /// still consumed, but no wake is reported.
    #[test]
    fn notify_without_a_consumer_reports_no_wake() {
        let d = MergeDemand::new();
        d.set_awaited(2);
        assert!(!d.notify_delivered(2));
        assert_eq!(d.awaited(), None);
        assert_eq!(d.snapshot().wakes_delivered, 0);
    }

    /// The unpark token: a notify that lands after `set_awaited` but before the
    /// consumer parks makes the park return at once. A lost token fails the 5 s
    /// watchdog (the park itself is 10 s). Between the notify and its park the
    /// consumer waits on an atomic with `yield_now`, not on a channel: a
    /// blocking `recv` parks the thread itself and would consume the token.
    #[test]
    fn notify_before_park_returns_immediately() {
        use std::sync::atomic::AtomicBool;
        let d = Arc::new(MergeDemand::new());
        let d2 = Arc::clone(&d);
        let notified = Arc::new(AtomicBool::new(false));
        let notified2 = Arc::clone(&notified);
        let (registered_tx, registered_rx) = channel::<()>();
        let (done_tx, done_rx) = channel::<()>();
        let h = std::thread::spawn(move || {
            d2.register_consumer();
            d2.set_awaited(1);
            registered_tx.send(()).unwrap();
            while !notified2.load(std::sync::atomic::Ordering::Acquire) {
                std::thread::yield_now();
            }
            std::thread::park_timeout(Duration::from_secs(10));
            done_tx.send(()).unwrap();
        });
        registered_rx.recv().unwrap();
        assert!(d.notify_delivered(1));
        notified.store(true, std::sync::atomic::Ordering::Release);
        done_rx
            .recv_timeout(Duration::from_secs(5))
            .expect("the park must return on the pending token");
        h.join().unwrap();
    }

    /// `await_slot` on a slot that already has a block (or EOF) clears `awaited`
    /// and returns `true`; on an empty, open slot it leaves `awaited` set. Each
    /// call is a registration; neither is a park until the caller records one.
    #[test]
    fn await_slot_rechecks_after_registering() {
        let d = MergeDemand::new();
        let slot = SortMergeSlot::new(
            4,
            std::io::BufReader::new(tempfile::tempfile().unwrap()),
            crate::codec::SpillCodec::Bgzf,
        );
        assert!(!d.await_slot(&slot), "empty, not EOF → park");
        assert_eq!(d.awaited(), Some(4));
        slot.decompressed.lock().unwrap().push_back(vec![1, 2, 3]);
        assert!(d.await_slot(&slot), "a block is there → do not park");
        assert_eq!(d.awaited(), None);
        let s = d.snapshot();
        assert_eq!((s.registrations, s.parking_registrations), (2, 0));
        d.stats().record_park();
        assert_eq!(d.snapshot().parking_registrations, 1);
    }

    /// `file_id == u32::MAX` must round-trip (the +1 encoding cannot overflow u64).
    #[test]
    fn max_file_id_round_trips() {
        let d = MergeDemand::new();
        d.set_awaited(u32::MAX);
        assert_eq!(d.awaited(), Some(u32::MAX));
        d.clear_awaited();
        assert_eq!(d.awaited(), None);
    }

    /// A different consumer thread replaces the stored handle.
    #[test]
    fn register_consumer_follows_the_current_thread() {
        let d = Arc::new(MergeDemand::new());
        d.register_consumer();
        let d2 = Arc::clone(&d);
        let other = std::thread::spawn(move || {
            d2.register_consumer();
            std::thread::current().id()
        })
        .join()
        .unwrap();
        assert_eq!(d.consumer_id_for_test(), Some(other));
    }

    #[test]
    fn starved_signal_is_shared_and_toggles() {
        let d = MergeDemand::new();
        let sig = d.starved_signal();
        assert!(!sig.load(std::sync::atomic::Ordering::Relaxed));
        d.set_starved(true);
        assert!(sig.load(std::sync::atomic::Ordering::Relaxed));
        d.set_starved(false);
        assert!(!sig.load(std::sync::atomic::Ordering::Relaxed));
        assert!(Arc::ptr_eq(&sig, &d.starved_signal()));
    }

    #[test]
    fn stall_buckets_and_log_lines() {
        let d = MergeDemand::new();
        d.stats().record_stall(AwaitedSlotState::Starved);
        d.stats().record_stall(AwaitedSlotState::Decompressing);
        d.stats().record_stall(AwaitedSlotState::Decompressing);
        d.stats().record_stall(AwaitedSlotState::Starved);
        d.stats().record_stall_ns(2_000_000_000);
        d.stats().record_partial_flush();
        assert!(d.stats().record_pool_request(MergePoolRequest::Woken, 3));
        assert!(!d.stats().record_pool_request(MergePoolRequest::AllAwake, 0));
        // Woken with nobody on the event-count: a timer-parked worker (or
        // worker 0 at one thread) was asleep, so it counts as a sleeper.
        assert!(d.stats().record_pool_request(MergePoolRequest::Woken, 0));
        assert!(!d.stats().record_pool_request(MergePoolRequest::Unavailable, 0));
        d.stats().record_stall_with_sleeper_ns(500_000_000);
        let s = d.snapshot();
        assert_eq!((s.stall_episodes, s.awaited_starved, s.awaited_decompressing), (4, 2, 2));
        assert_eq!(
            (s.pool_woken, s.pool_all_awake, s.pool_pending, s.pool_unavailable),
            (2, 1, 0, 1)
        );
        let lines = s.log_lines();
        assert_eq!(lines.len(), 4, "{lines:?}");
        assert!(lines[0].starts_with("Merge demand: 4 stall episodes (2.0 s exact)"), "{lines:?}");
        assert_eq!(lines[1], "Awaited slot at stall: starved 50% / decompressing 50%");
        assert_eq!(
            lines[2],
            "Pool at stall: requests woken 2 / all-awake 1 / pending 0 / unavailable 1; a \
             worker was asleep at 50.0% of requests (0.5 s)"
        );
        assert_eq!(lines[3], "Merge output: 1 partial flushes on stall");
        d.stats().mark_inline_decompress();
        let lines = d.snapshot().log_lines();
        assert_eq!(
            lines[1],
            "Awaited slot at stall: starved 50% / decompressing 50% \
             (inline decompress: stalls are not classified, every one reports starved)"
        );
    }
}

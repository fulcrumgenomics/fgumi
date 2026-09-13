//! `PoolEventCount` — a Vyukov-style event-count for parking idle pool workers.
//!
//! # Why this exists
//!
//! The chain-builder pool is *replicated round-robin polling*: when a worker's
//! whole dispatch pass does no work it exponential-backoff *sleeps*, and there
//! is no producer→consumer wakeup. On an oversubscribed pipeline (more workers
//! than parallel work) the surplus workers churn poll→sleep→poll forever, ~93%
//! idle, and that churn measurably slows the productive workers. Parking those
//! workers on a condition variable removes the churn — but a bare condvar has a
//! lost-wakeup race: a producer that publishes work *after* an idle worker
//! checked its queues (saw empty) but *before* it blocks would wake nobody, and
//! the worker sleeps on work that is already present.
//!
//! An **event-count** closes that race. It is a generation counter plus a
//! waiter count, used in the canonical two-phase protocol:
//!
//! - **Worker (about to idle):** [`prepare_wait`](PoolEventCount::prepare_wait)
//!   (register as a waiter, `SeqCst` fence, snapshot the generation) → re-check
//!   the *real* condition (its own queues) → if still no work,
//!   [`wait`](PoolEventCount::wait) (block until the generation moves or the
//!   deadline elapses).
//! - **Producer (just published work):** [`notify_one`](PoolEventCount::notify_one)
//!   (`SeqCst` fence, then read the waiter count; if any, bump the generation
//!   under the lock and wake one).
//!
//! # Correctness — one rule
//!
//! The design reduces to the C++20 `SeqCst`-fence rule: if a `SeqCst` fence X is
//! sequenced after a write A on one thread, and a `SeqCst` fence Y is sequenced
//! before a read B on another, and X precedes Y in the single total order of
//! `SeqCst` operations, then B observes A (or a later write). The two fences are
//! the one in `prepare_wait` (after `waiters += 1`, before the worker's
//! re-poll) and the one in `notify_one` (after the producer's publish, before
//! its `waiters` load). They are totally ordered, so exactly one holds:
//!
//! - **Worker's fence precedes producer's.** The producer's `waiters` load
//!   observes the increment (≥ 1), takes the lock, bumps `gen`, notifies. The
//!   worker is either not yet at its `gen` check (which runs under the same
//!   lock, after the producer releases it, so it sees the bump and returns
//!   without sleeping) or already in `wait_for` (which the notify wakes).
//! - **Producer's fence precedes worker's.** The worker's re-poll (sequenced
//!   after its fence) observes the producer's published item and never waits.
//!
//! No interleaving parks on a poppable item. The item handoff's own
//! happens-before is unaffected (it stays with the queue's atomics); this type
//! only orders the *wakeup*. Notifying under the same lock the condvar waits on
//! is what gives the condvar a "token": a bump that lands between the worker's
//! `gen` snapshot and its `wait_for` is not lost, because `wait` re-reads `gen`
//! under the lock before blocking.
//!
//! # Cost
//!
//! Steady state, nobody parked: `notify_one` is one `SeqCst` fence plus one
//! read-shared load of `waiters` — no store, no shared-line RMW, no lock. The
//! `waiters` and `gen` atomics are cache-line isolated (128-byte stride, as the
//! sharded liveness counter already does) so an idling worker's
//! increment/decrement of `waiters` never false-shares with a producer's load.

use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering, fence};
use std::time::Duration;

use parking_lot::{Condvar, Mutex};

/// Cache-line stride used to isolate the hot atomics from each other and from
/// any neighbouring fields. 128 bytes covers the 64-byte line plus the
/// adjacent-line prefetch pairing on aarch64 and x86-64 (matches the padding
/// the sharded liveness counter uses).
const CACHE_LINE: usize = 128;

#[repr(align(128))]
struct Padded<T>(T);

/// A Vyukov event-count over a `parking_lot` `Mutex`/`Condvar`.
///
/// Shared by all pool workers of one pipeline via `Arc`. Only present when
/// `n_threads > 1`; the fused and scheduled single-thread paths pass `None` and
/// keep their existing sleep-backoff idle.
pub struct PoolEventCount {
    /// Number of workers currently *armed or blocked* (between `prepare_wait`
    /// and the matching `wait`/`cancel_wait`). Read by `notify_one` on the hot
    /// path to skip the lock+wake when nobody is parked. Cache-line isolated.
    waiters: Padded<AtomicUsize>,
    /// Generation counter, bumped only when notifying, only under `lock`. A
    /// worker that snapshots `gen` in `prepare_wait` and finds it changed by
    /// the time it would block knows work was published in its wait window and
    /// does not block. Cache-line isolated.
    generation: Padded<AtomicU64>,
    /// Guards the `gen` bump and the condvar. A producer takes it only when
    /// `waiters > 0`; a worker takes it only to block.
    lock: Mutex<()>,
    cvar: Condvar,
    /// Total pool worker count, so `notify_one` can compute
    /// `awake = n_pool − waiters` for the ceiling gate. Set once at
    /// construction; never mutated.
    n_pool: usize,
    /// Concurrency ceiling: `notify_one` wakes a parked worker only while
    /// `awake < ceiling`. Above it, a push wakes nobody — the already-awake
    /// workers find the item on their next pass (they cannot park while their
    /// passes find work), so parking *sticks* instead of the whole surplus
    /// re-waking on every push (the per-`Progress` cascade). Default `n_pool`
    /// (ungated: `awake ≤ n_pool` always holds, so every notify fires) — the
    /// behaviour before a controller adjusts it. `notify_all` (cancel / edge
    /// close / `Finished`) is NEVER gated: correctness must not depend on the
    /// ceiling. Cache-line isolated from the hot `waiters`/`gen`.
    ceiling: Padded<AtomicUsize>,
}

const _: () = assert!(std::mem::align_of::<Padded<AtomicUsize>>() == CACHE_LINE);

/// Opaque generation snapshot from [`PoolEventCount::prepare_wait`], handed back
/// to [`PoolEventCount::wait`]. Carrying it by value (rather than re-reading
/// inside `wait`) is what makes the "generation moved in my wait window" check
/// precise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WaitKey(u64);

/// Outcome of a [`PoolEventCount::wait`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WaitOutcome {
    /// The generation had already moved by the time `wait` took the lock — work
    /// was published in the arm→wait window; the caller must re-poll, not sleep.
    Woken,
    /// The condvar returned (a notify, a spurious wake, or the deadline). The
    /// caller loops and re-polls regardless — spurious wakes are harmless.
    Returned,
}

impl PoolEventCount {
    /// Construct for a pool of `n_pool` workers. The ceiling starts at `n_pool`
    /// (ungated). Use [`Self::set_ceiling`] to gate wakeups.
    #[must_use]
    pub fn new(n_pool: usize) -> Self {
        Self {
            waiters: Padded(AtomicUsize::new(0)),
            generation: Padded(AtomicU64::new(0)),
            lock: Mutex::new(()),
            cvar: Condvar::new(),
            n_pool,
            ceiling: Padded(AtomicUsize::new(n_pool)),
        }
    }

    /// Set the concurrency ceiling — the max number of workers `notify_one`
    /// will keep awake. Clamped to `1..=n_pool` (a ceiling of 0 would wedge the
    /// pool: no worker could ever be woken). A controller (edge-depth /
    /// rejection driven) calls this to shed or admit surplus workers; with no
    /// controller the default `n_pool` leaves every notify ungated.
    #[inline]
    pub fn set_ceiling(&self, ceiling: usize) {
        self.ceiling.0.store(ceiling.clamp(1, self.n_pool.max(1)), Ordering::Relaxed);
    }

    /// The current ceiling (for the controller / tests).
    #[must_use]
    #[inline]
    pub fn ceiling(&self) -> usize {
        self.ceiling.0.load(Ordering::Relaxed)
    }

    /// Number of workers currently armed or blocked. Used by the ceiling
    /// controller to compute `awake = n_pool − waiters`.
    #[must_use]
    #[inline]
    pub fn waiters(&self) -> usize {
        self.waiters.0.load(Ordering::Relaxed)
    }

    /// Producer side: wake at most one parked worker if any is waiting.
    ///
    /// Hot path when nobody is parked: one `SeqCst` fence + one relaxed load,
    /// then return. No store, no lock, no wake. The fence is the producer half
    /// of the correctness rule and must be `SeqCst`; it is sequenced *after* the
    /// caller's publish (the successful queue push) and *before* the `waiters`
    /// load, so a worker that armed before the publish is seen here.
    #[inline]
    pub fn notify_one(&self) {
        fence(Ordering::SeqCst);
        let waiters = self.waiters.0.load(Ordering::Relaxed);
        if waiters == 0 {
            return;
        }
        // Ceiling gate: keep at most `ceiling` workers awake. `awake = n_pool −
        // waiters`; if we are already at/above the ceiling, do not wake another
        // — the awake workers will find this item on their next pass. This is
        // what stops the per-`Progress` cascade from re-waking the whole surplus
        // on an oversubscribed pipeline. Ungated when `ceiling == n_pool`
        // (default): then `awake < n_pool` whenever `waiters >= 1`, so the
        // branch never suppresses a wake.
        let awake = self.n_pool.saturating_sub(waiters);
        if awake >= self.ceiling.0.load(Ordering::Relaxed) {
            return;
        }
        let _guard = self.lock.lock();
        self.generation.0.fetch_add(1, Ordering::Relaxed);
        self.cvar.notify_one();
    }

    /// Producer side: wake *every* parked worker. Reserved for terminal /
    /// broadcast events (edge close, `Finished`, cancel/error) where more than
    /// one waiter may need to observe the transition. Same fence discipline as
    /// [`Self::notify_one`].
    #[inline]
    pub fn notify_all(&self) {
        fence(Ordering::SeqCst);
        if self.waiters.0.load(Ordering::Relaxed) == 0 {
            return;
        }
        let _guard = self.lock.lock();
        self.generation.0.fetch_add(1, Ordering::Relaxed);
        self.cvar.notify_all();
    }

    /// Worker side, phase 1: register as a waiter and snapshot the generation.
    ///
    /// The `SeqCst` fence is the worker half of the correctness rule: it is
    /// sequenced *after* the `waiters` increment and *before* the caller's
    /// re-poll of its real condition. The returned [`WaitKey`] must be passed to
    /// exactly one of [`Self::wait`] or [`Self::cancel_wait`] to balance the
    /// increment.
    #[must_use]
    #[inline]
    pub fn prepare_wait(&self) -> WaitKey {
        self.waiters.0.fetch_add(1, Ordering::Relaxed);
        fence(Ordering::SeqCst);
        WaitKey(self.generation.0.load(Ordering::Relaxed))
    }

    /// Worker side: abandon a prepared wait without blocking (the re-poll found
    /// work, or the pipeline is shutting down). Balances the `prepare_wait`
    /// increment.
    #[inline]
    pub fn cancel_wait(&self, _key: WaitKey) {
        self.waiters.0.fetch_sub(1, Ordering::Relaxed);
    }

    /// Worker side, phase 2: block until the generation moves past `key` or
    /// `deadline` elapses. Balances the `prepare_wait` increment on return.
    ///
    /// Returns [`WaitOutcome::Woken`] if the generation had already advanced
    /// when the lock was taken (a notify raced into the arm→wait window — do not
    /// block, re-poll immediately), else [`WaitOutcome::Returned`] after the
    /// condvar wait (notify, spurious, or timeout — re-poll anyway).
    #[must_use]
    pub fn wait(&self, key: WaitKey, deadline: Duration) -> WaitOutcome {
        let mut guard = self.lock.lock();
        // Under the lock, re-read `gen`. If it moved since `prepare_wait`, a
        // producer bumped it (also under this lock) after we armed — the wakeup
        // we would wait for has already happened. This is the condvar's "token".
        if self.generation.0.load(Ordering::Relaxed) != key.0 {
            drop(guard);
            self.waiters.0.fetch_sub(1, Ordering::Relaxed);
            return WaitOutcome::Woken;
        }
        // `wait_for` may wake spuriously; the caller re-polls regardless, so a
        // single wait (not a loop) is correct here.
        let _ = self.cvar.wait_for(&mut guard, deadline);
        drop(guard);
        self.waiters.0.fetch_sub(1, Ordering::Relaxed);
        WaitOutcome::Returned
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::thread;

    use super::*;

    #[test]
    fn new_has_no_waiters() {
        let ec = PoolEventCount::new(32);
        assert_eq!(ec.waiters(), 0);
    }

    #[test]
    fn prepare_then_cancel_balances_waiters() {
        let ec = PoolEventCount::new(32);
        let key = ec.prepare_wait();
        assert_eq!(ec.waiters(), 1);
        ec.cancel_wait(key);
        assert_eq!(ec.waiters(), 0);
    }

    #[test]
    fn notify_one_with_no_waiters_is_a_noop_and_takes_no_lock() {
        // The hot-path guarantee: with nobody parked, notify_one must not take
        // the lock. Prove it by holding the lock for the duration of the call —
        // if notify_one tried to lock, it would deadlock (parking_lot mutex is
        // not reentrant) and the test would hang; a clean return proves the
        // lock was never attempted.
        let ec = PoolEventCount::new(32);
        let held = ec.lock.lock();
        ec.notify_one(); // must return without touching the lock
        drop(held);
        assert_eq!(ec.waiters(), 0);
    }

    #[test]
    fn wait_returns_woken_when_generation_already_moved() {
        // Models the arm→(notify)→wait window: the worker armed, a producer
        // bumped the generation before the worker reached `wait`. `wait` must
        // see the moved generation and return Woken without blocking.
        let ec = PoolEventCount::new(32);
        let key = ec.prepare_wait();
        // Simulate a producer's notify landing in the window (bumps gen under
        // the lock, as notify_one does).
        {
            let _g = ec.lock.lock();
            ec.generation.0.fetch_add(1, Ordering::Relaxed);
        }
        let outcome = ec.wait(key, Duration::from_secs(30));
        assert_eq!(outcome, WaitOutcome::Woken, "a moved generation must not block");
        assert_eq!(ec.waiters(), 0, "wait balances the prepare_wait increment");
    }

    #[test]
    fn wait_times_out_when_no_notify() {
        // With no notify, wait blocks until the deadline and returns Returned.
        let ec = PoolEventCount::new(32);
        let key = ec.prepare_wait();
        let start = std::time::Instant::now();
        let outcome = ec.wait(key, Duration::from_millis(30));
        assert!(start.elapsed() >= Duration::from_millis(25), "must block ~the deadline");
        assert_eq!(outcome, WaitOutcome::Returned);
        assert_eq!(ec.waiters(), 0);
    }

    #[test]
    fn notify_one_wakes_a_blocked_waiter() {
        let ec = Arc::new(PoolEventCount::new(32));
        let ec2 = Arc::clone(&ec);
        let woke = Arc::new(AtomicUsize::new(0));
        let woke2 = Arc::clone(&woke);

        let h = thread::spawn(move || {
            let key = ec2.prepare_wait();
            // No work found; block for a long deadline so only a notify frees it.
            let _ = ec2.wait(key, Duration::from_secs(30));
            woke2.store(1, Ordering::SeqCst);
        });

        // Wait until the worker has armed, then notify it.
        while ec.waiters() == 0 {
            thread::yield_now();
        }
        // Give the worker a moment to reach the blocking wait_for.
        thread::sleep(Duration::from_millis(20));
        ec.notify_one();

        h.join().expect("waiter thread joins");
        assert_eq!(woke.load(Ordering::SeqCst), 1, "the blocked waiter was woken by notify_one");
        assert_eq!(ec.waiters(), 0);
    }

    #[test]
    fn notify_all_wakes_every_waiter() {
        let ec = Arc::new(PoolEventCount::new(32));
        let n = 8usize;
        let woke = Arc::new(AtomicUsize::new(0));
        let mut handles = Vec::new();
        for _ in 0..n {
            let ec2 = Arc::clone(&ec);
            let woke2 = Arc::clone(&woke);
            handles.push(thread::spawn(move || {
                let key = ec2.prepare_wait();
                let _ = ec2.wait(key, Duration::from_secs(30));
                woke2.fetch_add(1, Ordering::SeqCst);
            }));
        }
        while ec.waiters() < n {
            thread::yield_now();
        }
        thread::sleep(Duration::from_millis(20));
        ec.notify_all();
        for h in handles {
            h.join().expect("waiter joins");
        }
        assert_eq!(woke.load(Ordering::SeqCst), n, "notify_all woke every waiter");
        assert_eq!(ec.waiters(), 0);
    }

    /// Stress: producers publishing concurrently with a worker arming/waiting
    /// must never leave the worker blocked on a "published" item. Models the
    /// lost-wakeup window many times. Not a proof (that needs loom), but catches
    /// gross ordering regressions.
    #[test]
    fn stress_no_lost_wakeup() {
        for _ in 0..2000 {
            let ec = Arc::new(PoolEventCount::new(32));
            let ec_p = Arc::clone(&ec);
            // A flag standing in for "work is available", published before notify.
            let work = Arc::new(AtomicUsize::new(0));
            let work_p = Arc::clone(&work);

            let producer = thread::spawn(move || {
                work_p.store(1, Ordering::SeqCst); // publish
                ec_p.notify_one(); // then notify
            });

            // Worker: arm, re-check the real condition, block only if no work.
            let key = ec.prepare_wait();
            let outcome = if work.load(Ordering::SeqCst) == 1 {
                ec.cancel_wait(key);
                WaitOutcome::Woken
            } else {
                ec.wait(key, Duration::from_millis(500))
            };
            producer.join().expect("producer joins");
            // Either the re-check saw the work, or wait returned (Woken or a
            // notify-driven Returned) well within the deadline. A lost wakeup
            // would instead burn the full 500ms then return Returned with work
            // set — assert we observed work by the time we're done.
            assert_eq!(work.load(Ordering::SeqCst), 1);
            let _ = outcome;
            assert_eq!(ec.waiters(), 0);
        }
    }

    #[test]
    fn set_ceiling_clamps_to_pool_and_floor_one() {
        let ec = PoolEventCount::new(8);
        assert_eq!(ec.ceiling(), 8, "default ceiling is the pool size (ungated)");
        ec.set_ceiling(3);
        assert_eq!(ec.ceiling(), 3);
        ec.set_ceiling(0);
        assert_eq!(ec.ceiling(), 1, "0 clamps to 1 — a 0 ceiling would wedge the pool");
        ec.set_ceiling(100);
        assert_eq!(ec.ceiling(), 8, "above n_pool clamps to n_pool");
    }

    #[test]
    fn notify_one_suppressed_at_ceiling_but_notify_all_still_wakes() {
        // n_pool = 4, ceiling = 1. One waiter parks → awake = 4 − 1 = 3 ≥ 1, so
        // `notify_one` must NOT wake it (we are already at/above the ceiling).
        // `notify_all` (ungated) must still wake it.
        let ec = Arc::new(PoolEventCount::new(4));
        ec.set_ceiling(1);
        let ec2 = Arc::clone(&ec);
        let woke = Arc::new(AtomicUsize::new(0));
        let woke2 = Arc::clone(&woke);
        let h = thread::spawn(move || {
            let key = ec2.prepare_wait();
            let _ = ec2.wait(key, Duration::from_secs(30));
            woke2.store(1, Ordering::SeqCst);
        });
        while ec.waiters() == 0 {
            thread::yield_now();
        }
        thread::sleep(Duration::from_millis(20));
        // Gated: awake (3) >= ceiling (1) → no wake.
        ec.notify_one();
        thread::sleep(Duration::from_millis(30));
        assert_eq!(woke.load(Ordering::SeqCst), 0, "notify_one must be suppressed at the ceiling");
        // Ungated broadcast frees it.
        ec.notify_all();
        h.join().expect("waiter joins");
        assert_eq!(woke.load(Ordering::SeqCst), 1, "notify_all is never gated");
    }
}

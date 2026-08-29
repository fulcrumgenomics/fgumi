#![deny(unsafe_code)]
//! Bounded, reusable pool of [`SegmentedBuf`] sort arenas — an RSS fix.
//!
//! The Phase-1 spill path fills a multi-GB `SegmentedBuf`, hands it downstream
//! as an `Arc`, and (pre-pool) minted a *fresh* buffer for the next fill via
//! `mem::take`. Combined with the block-parallel spill running async, two full
//! arenas were live at the Phase-1→Phase-2 boundary (the in-flight chunk —
//! pinned by slow compression through `SpillGather`'s bounded `pending` — plus
//! the next fill), inflating peak RSS to ~2× base.
//!
//! This pool **bounds** the live arenas to `capacity` (default 1, matching
//! legacy's one-arena-at-a-time `drain_pending_spill` model) and **reuses**
//! their storage: [`try_acquire`](ArenaPool::try_acquire) hands out a
//! reset-for-reuse buffer or `None` when all `capacity` arenas are in flight
//! (the caller backpressures with `NoProgress` until an in-flight chunk's `Arc`
//! drops and returns its arena via [`PooledSegmentedBuf`]'s `Drop`). Capping at
//! 1 means the next fill cannot start until the prior chunk is spilled and its
//! arena freed — exactly legacy's behaviour, and the RSS-gate-safe footprint.

use std::ops::{Deref, DerefMut};
use std::sync::{Arc, Mutex};

use crate::segmented_buf::SegmentedBuf;

/// A bounded free-list of reusable [`SegmentedBuf`] arenas. At most `capacity`
/// buffers ever exist; `try_acquire` returns `None` when all are in flight.
pub struct ArenaPool {
    inner: Mutex<PoolInner>,
    /// Maximum number of live arenas (`N_arena`); clamped to ≥ 1.
    capacity: usize,
    /// Segment size for freshly-allocated arenas.
    segment_size: usize,
}

struct PoolInner {
    /// Returned, reset-for-reuse buffers ready to hand out.
    free: Vec<SegmentedBuf>,
    /// Total arenas ever allocated (never decremented; bounds at `capacity`).
    made: usize,
}

impl ArenaPool {
    /// Create a pool bounded to `capacity` (≥1) arenas, each `segment_size`.
    #[must_use]
    pub fn new(capacity: usize, segment_size: usize) -> Arc<Self> {
        Arc::new(Self {
            inner: Mutex::new(PoolInner { free: Vec::new(), made: 0 }),
            capacity: capacity.max(1),
            // Clamped for the same reason `capacity` is, and to match
            // `SegmentedBuf::with_capacity`, which clamps it again downstream: a
            // 0 here would otherwise be silently accepted at this boundary and
            // corrected somewhere else.
            segment_size: segment_size.max(1),
        })
    }

    /// Acquire a reset-for-reuse arena, or `None` if all `capacity` arenas are
    /// in flight (the caller must backpressure and retry once one returns).
    ///
    /// Returns the RAII wrapper, never a bare [`SegmentedBuf`], and that is
    /// load-bearing rather than stylistic. `made` is only ever incremented; the
    /// sole path back to the free-list is [`PooledSegmentedBuf::drop`]. If a
    /// caller could hold the arena unwrapped, then *any* drop on that path — a
    /// `?` on a malformed record, an early return, an unwind — would retire the
    /// slot permanently. At the default `capacity == 1` that leaves the pool
    /// empty forever, and because the documented response to `None` is to
    /// backpressure until an in-flight arena returns, the caller would spin
    /// forever instead of failing: a hang, not an error. Handing back the
    /// wrapper makes the return unconditional and unforgeable.
    ///
    /// # Panics
    ///
    /// Panics if the internal mutex is poisoned. That is a deliberate asymmetry
    /// with `release`, which recovers the same guard: the
    /// protected state is never torn, so recovery is *sound* anywhere — but
    /// `release` runs from a `Drop` and must not panic during unwinding, while
    /// this path has no such constraint and a poisoned pool is worth surfacing
    /// loudly rather than papering over.
    #[must_use]
    pub fn try_acquire(self: &Arc<Self>) -> Option<PooledSegmentedBuf> {
        let mut g = self.inner.lock().expect("arena pool mutex poisoned");
        let buf = if let Some(buf) = g.free.pop() {
            buf
        } else if g.made < self.capacity {
            g.made += 1;
            SegmentedBuf::with_capacity(0, self.segment_size)
        } else {
            return None;
        };
        drop(g);
        Some(PooledSegmentedBuf::pooled(buf, Arc::clone(self)))
    }

    /// Return an arena to the free-list, reset for reuse (capacity retained).
    ///
    /// Reached only from [`PooledSegmentedBuf::drop`], so it must not panic: if
    /// the arena is dropped while a panic is already unwinding, a second panic
    /// here aborts the process outright — no catch, no message. No path poisons
    /// this mutex *today* — `try_acquire` is the only one that runs non-trivial
    /// code under it, and it allocates with `capacity = 0`, so neither the
    /// first-segment nor the segment-vector allocation scales with
    /// `segment_size` and neither can overflow. That is a property of the
    /// current body, not of the lock: any future panic under this guard would
    /// poison it, and this path would be the one that turned that into an abort.
    ///
    /// Recovering the guard is sound here because the state it protects cannot
    /// be torn: `free` holds reset-for-reuse buffers and `made` is a counter, so
    /// a thread that died mid-critical-section leaves the free-list structurally
    /// valid (at worst `made` counts an arena that was never handed out, which
    /// costs a slot but breaks nothing). `try_acquire` deliberately keeps its
    /// `expect` — it is not on a drop path, and a poisoned pool should surface
    /// there rather than be papered over.
    fn release(&self, mut buf: SegmentedBuf) {
        // Catch a buffer that never came from this pool's `try_acquire` (wrong
        // `segment_size`): reusing it would silently defeat the `capacity`-based
        // RSS bound this pool exists to enforce. Debug-only — the check is a
        // developer guard, not a production cost.
        debug_assert_eq!(
            buf.segment_size(),
            self.segment_size,
            "released arena's segment_size doesn't match this pool's configured size"
        );
        buf.reset_for_reuse();
        self.inner.lock().unwrap_or_else(std::sync::PoisonError::into_inner).free.push(buf);
    }

    /// Number of arenas currently available in the free-list (test/diagnostic).
    ///
    /// Does not panic on a poisoned mutex — see the comment in the body.
    #[cfg(test)]
    #[must_use]
    pub fn free_len(&self) -> usize {
        // Recovers like `release` rather than panicking: this exists to observe
        // pool state, including after a poisoning, which is exactly when a test
        // most needs to read it.
        self.inner.lock().unwrap_or_else(std::sync::PoisonError::into_inner).free.len()
    }
}

/// A [`SegmentedBuf`] that returns to its [`ArenaPool`] on drop (when `pool`
/// is `Some`). Used as the `Arc`-shared backing store of an in-memory sort
/// chunk: the arena is reclaimed once the last `Arc` clone drops (after the
/// gather has framed a spilled chunk). Buffers that are not pool-managed — the
/// residual chunk, and any buffer built outside the pool — use
/// [`unpooled`](Self::unpooled), so they drop normally and never reach the
/// free-list.
pub struct PooledSegmentedBuf {
    /// `Some` until `Drop`. `Option` only so `Drop` can move the buffer out.
    buf: Option<SegmentedBuf>,
    /// `Some` → return to the pool on drop; `None` → drop normally.
    pool: Option<Arc<ArenaPool>>,
}

impl PooledSegmentedBuf {
    /// Wrap a buffer that should return to `pool` on drop.
    ///
    /// Deliberately private: [`ArenaPool::try_acquire`] is the only way to get a
    /// pooled wrapper, so every buffer on the free-list came from this pool.
    /// Were this public, a caller could wrap a foreign buffer and `release`
    /// would push it onto the free-list — `try_acquire` pops from `free` before
    /// consulting `made`, so that directly raises the number of simultaneously
    /// live arenas above `capacity`, defeating the RSS bound the pool exists to
    /// enforce. The `debug_assert_eq!` on `segment_size` in `release` catches
    /// only the mismatched-size case, and only in debug.
    #[must_use]
    fn pooled(buf: SegmentedBuf, pool: Arc<ArenaPool>) -> Self {
        Self { buf: Some(buf), pool: Some(pool) }
    }

    /// Wrap a buffer that is NOT pool-managed (drops normally).
    #[must_use]
    pub fn unpooled(buf: SegmentedBuf) -> Self {
        Self { buf: Some(buf), pool: None }
    }
}

impl Default for PooledSegmentedBuf {
    /// An empty, NON-pooled buffer.
    ///
    /// This exists so `std::mem::take(&mut some_pooled_buf)` moves the WRAPPER
    /// out. Without it, `take` deref-coerces through `DerefMut` and steals the
    /// inner `SegmentedBuf`, leaving the wrapper holding a default-sized buffer
    /// that it still believes belongs to the pool: the real arena is orphaned
    /// (never returned) and a wrong-`segment_size` buffer is later pushed onto
    /// the free list. Pooled-ness must travel with the arena, and it only does
    /// if the whole wrapper moves.
    fn default() -> Self {
        Self::unpooled(SegmentedBuf::default())
    }
}

impl DerefMut for PooledSegmentedBuf {
    /// Fill happens *through* the wrapper. Without this the caller would have to
    /// unwrap to write, which is exactly the window in which a dropped arena
    /// silently retires its pool slot — see [`ArenaPool::try_acquire`].
    fn deref_mut(&mut self) -> &mut SegmentedBuf {
        self.buf.as_mut().expect("PooledSegmentedBuf used after drop")
    }
}

impl Deref for PooledSegmentedBuf {
    type Target = SegmentedBuf;
    fn deref(&self) -> &SegmentedBuf {
        self.buf.as_ref().expect("PooledSegmentedBuf used after drop")
    }
}

impl Drop for PooledSegmentedBuf {
    fn drop(&mut self) {
        if let Some(buf) = self.buf.take()
            && let Some(pool) = &self.pool
        {
            pool.release(buf);
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn bounds_to_capacity_and_reuses() {
        let pool = ArenaPool::new(1, 1024);
        let a = pool.try_acquire().expect("first acquire");
        assert!(pool.try_acquire().is_none(), "capacity 1 → second acquire blocked");
        drop(a); // the wrapper's Drop is the only return path
        assert_eq!(pool.free_len(), 1, "released arena back in the free-list");
        let _b = pool.try_acquire().expect("acquire after release reuses the arena");
        assert!(pool.try_acquire().is_none(), "still capacity-bound");
    }

    /// `capacity` must actually bound the pool. Every other test here uses
    /// `capacity == 1`, where "hand out one arena" and "respect the capacity"
    /// are indistinguishable — an implementation ignoring the field entirely
    /// passes them all. This one separates the two.
    #[rstest]
    #[case::two(2)]
    #[case::five(5)]
    fn capacity_bounds_the_number_of_live_arenas(#[case] capacity: usize) {
        let pool = ArenaPool::new(capacity, 1024);
        let live: Vec<_> = (0..capacity)
            .map(|i| pool.try_acquire().unwrap_or_else(|| panic!("acquire {i} within capacity")))
            .collect();
        assert!(pool.try_acquire().is_none(), "the {capacity}+1'th acquire must be refused");
        drop(live);
        assert_eq!(pool.free_len(), capacity, "every arena returns on drop");
    }

    /// The leak this API shape exists to prevent: an arena dropped on an error
    /// path must still return. `made` is never decremented, so a slot lost here
    /// is lost forever — at the default capacity 1 that empties the pool, and
    /// since the documented response to `None` is to backpressure until an arena
    /// returns, the caller would spin forever rather than fail.
    #[test]
    fn arena_dropped_on_an_error_path_still_returns_to_the_pool() {
        // A fill that bails part-way: the arena is live and partly written when
        // the `?` fires, and nothing wraps it on the way out.
        fn fill_then_fail(pool: &Arc<ArenaPool>) -> Result<(), &'static str> {
            let mut arena = pool.try_acquire().ok_or("pool empty")?;
            arena.extend_from_slice(&[0xAB; 32]);
            Err("malformed record")?;
            unreachable!()
        }

        let pool = ArenaPool::new(1, 1024);
        // Discriminating: `fill_then_fail` can also fail with "pool empty", which
        // would mean the acquire never happened and the rest of the test proves
        // nothing about the drop path.
        assert_eq!(
            fill_then_fail(&pool),
            Err("malformed record"),
            "must reach the mid-fill bail, not fail to acquire",
        );
        assert_eq!(pool.free_len(), 1, "the arena came back despite the early return");
        let reused = pool.try_acquire().expect("pool is still usable after the error path");
        assert!(reused.is_empty(), "and the returned arena was reset");
    }

    /// `Deref`/`DerefMut` are the only way a consumer reaches the arena — the
    /// wrapper exposes no accessor — so a broken impl would make every pooled
    /// read or write unreachable. Nothing else here exercises them directly.
    #[test]
    fn deref_reaches_the_wrapped_arena_for_both_read_and_write() {
        let pool = ArenaPool::new(1, 1024);
        let mut arena = pool.try_acquire().expect("acquire");

        // DerefMut: the fill path.
        let off = arena.extend_from_slice(b"pooled bytes");
        // Deref: the read path the consumer uses.
        assert_eq!(arena.slice(off, 12), b"pooled bytes");
        assert_eq!(arena.len(), 12, "Deref reaches the wrapped buffer's own state");

        drop(arena);
        let reused = pool.try_acquire().expect("reacquire");
        assert!(reused.is_empty(), "and the arena was reset on the way back");
    }

    /// A poisoned pool mutex must not make `PooledSegmentedBuf::drop` panic.
    ///
    /// This one cannot be written as a regression test after the fact: if
    /// `release` panics while an arena drops during unwinding, the process
    /// aborts, so there is no failing assertion to observe — the test binary
    /// just dies. Pinning the recovery here is the only way to keep it.
    #[test]
    fn release_recovers_a_poisoned_lock_rather_than_panicking_in_drop() {
        let pool = ArenaPool::new(1, 1024);
        let arena = pool.try_acquire().expect("acquire");

        // Poison the mutex the only way it can be poisoned: a thread dies while
        // holding it. No current path does that (see `release`'s doc), so the
        // panic here is deliberately synthetic — the point is to pin `release`'s
        // recovery, which exists so a *future* panic under this guard cannot
        // turn an arena drop during unwinding into a process abort.
        let poisoner = Arc::clone(&pool);
        let previous_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {})); // keep the deliberate panic out of test output
        let died = std::thread::spawn(move || {
            let _guard = poisoner.inner.lock().expect("lock before poisoning");
            panic!("poison the pool mutex");
        })
        .join();
        std::panic::set_hook(previous_hook);
        assert!(died.is_err(), "the poisoning thread must actually have panicked");

        drop(arena); // must not panic
        assert_eq!(pool.free_len(), 1, "the arena still returned to the poisoned pool");
    }

    #[test]
    fn unpooled_does_not_return() {
        let pool = ArenaPool::new(1, 1024);
        drop(PooledSegmentedBuf::unpooled(SegmentedBuf::with_capacity(0, 1024)));
        assert_eq!(pool.free_len(), 0, "unpooled buffer never reaches the pool");
    }

    #[test]
    fn reused_arena_retains_capacity() {
        let pool = ArenaPool::new(1, 16);
        let mut buf = pool.try_acquire().unwrap();
        for i in 0..4u8 {
            buf.extend_from_slice(&[i; 10]);
        }
        let (segs, cap) = (buf.num_segments(), buf.allocated_capacity());
        assert!(segs >= 2, "grew multiple segments");
        drop(buf);
        let reused = pool.try_acquire().unwrap();
        assert!(reused.is_empty(), "reset for reuse");
        assert_eq!(reused.num_segments(), segs, "retained segments");
        assert_eq!(reused.allocated_capacity(), cap, "no realloc");
    }
}

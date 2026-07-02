#![deny(unsafe_code)]
//! Bounded, reusable pool of [`SegmentedBuf`] sort arenas (lever-1 RSS fix).
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

use std::ops::Deref;
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
            segment_size,
        })
    }

    /// Acquire a reset-for-reuse arena, or `None` if all `capacity` arenas are
    /// in flight (the caller must backpressure and retry once one returns).
    ///
    /// # Panics
    ///
    /// Panics if the internal mutex is poisoned (another thread panicked while
    /// holding it, which indicates an unrecoverable state).
    #[must_use]
    pub fn try_acquire(&self) -> Option<SegmentedBuf> {
        let mut g = self.inner.lock().expect("arena pool mutex poisoned");
        if let Some(buf) = g.free.pop() {
            return Some(buf);
        }
        if g.made < self.capacity {
            g.made += 1;
            return Some(SegmentedBuf::with_capacity(0, self.segment_size));
        }
        None
    }

    /// Return an arena to the free-list, reset for reuse (capacity retained).
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
        self.inner.lock().expect("arena pool mutex poisoned").free.push(buf);
    }

    /// Number of arenas currently available in the free-list (test/diagnostic).
    ///
    /// # Panics
    ///
    /// Panics if the internal mutex is poisoned.
    #[cfg(test)]
    #[must_use]
    pub fn free_len(&self) -> usize {
        self.inner.lock().expect("arena pool mutex poisoned").free.len()
    }
}

/// A [`SegmentedBuf`] that returns to its [`ArenaPool`] on drop (when `pool`
/// is `Some`). Used as the `Arc`-shared backing store of an in-memory sort
/// chunk: the arena is reclaimed once the last `Arc` clone drops (after the
/// gather has framed a spilled chunk). The residual chunk and the
/// `from_owned_records`/`empty` test constructors use `pool == None` so they
/// drop normally and never return to the pool.
pub struct PooledSegmentedBuf {
    /// `Some` until `Drop`. `Option` only so `Drop` can move the buffer out.
    buf: Option<SegmentedBuf>,
    /// `Some` → return to the pool on drop; `None` → drop normally.
    pool: Option<Arc<ArenaPool>>,
}

impl PooledSegmentedBuf {
    /// Wrap a buffer that should return to `pool` on drop.
    #[must_use]
    pub fn pooled(buf: SegmentedBuf, pool: Arc<ArenaPool>) -> Self {
        Self { buf: Some(buf), pool: Some(pool) }
    }

    /// Wrap a buffer that is NOT pool-managed (drops normally).
    #[must_use]
    pub fn unpooled(buf: SegmentedBuf) -> Self {
        Self { buf: Some(buf), pool: None }
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
        if let Some(buf) = self.buf.take() {
            if let Some(pool) = &self.pool {
                pool.release(buf);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_to_capacity_and_reuses() {
        let pool = ArenaPool::new(1, 1024);
        let a = pool.try_acquire().expect("first acquire");
        assert!(pool.try_acquire().is_none(), "capacity 1 → second acquire blocked");
        // Return it via the Pooled wrapper's Drop.
        drop(PooledSegmentedBuf::pooled(a, Arc::clone(&pool)));
        assert_eq!(pool.free_len(), 1, "released arena back in the free-list");
        let _b = pool.try_acquire().expect("acquire after release reuses the arena");
        assert!(pool.try_acquire().is_none(), "still capacity-bound");
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
        pool.release(buf);
        let reused = pool.try_acquire().unwrap();
        assert!(reused.is_empty(), "reset for reuse");
        assert_eq!(reused.num_segments(), segs, "retained segments");
        assert_eq!(reused.allocated_capacity(), cap, "no realloc");
    }
}

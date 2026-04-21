//! Global memory accounting via atomic CAS.
//!
//! # Scope of the counter
//!
//! This tracker measures **bytes held in `StageQueue`s** at any given moment.
//! The byte count per item comes from [`crate::runall::engine::stage::Stage::output_memory_estimate`],
//! which each stage implements to return an estimate of the item's heap usage
//! (e.g. `vec.capacity()` for a `Vec<u8>` payload). The counter is increased
//! on a successful `StageQueue::push` and decreased on a successful
//! `StageQueue::pop`.
//!
//! What the counter does NOT track:
//! - Payloads held in sink / coalesce `ReorderBuffer`s after being popped.
//! - Items stashed in worker held-item slots after a failed push.
//! - Items currently being processed inside a stage's `process()` call
//!   (popped from the input queue but not yet pushed to the output queue).
//! - Stage-internal scratch buffers (decompression arenas, SAM parser state,
//!   ring buffers).
//! - Subprocess memory (e.g. the aligner's `bwa mem` child process).
//! - Allocator overhead, fragmentation, thread stacks.
//!
//! As a result, `MemoryTracker::peak()` is a lower bound on process RSS, not
//! an estimate of it. See the logs emitted by [`crate::runall::engine::driver`] for
//! the user-facing wording.

use std::sync::atomic::{AtomicUsize, Ordering};

/// Tracks total bytes across all pipeline queues.
///
/// Uses a lock-free compare-and-swap loop so concurrent producers and consumers
/// can update the counter without a mutex. See the module-level docs for the
/// scope of what this counter does and does not track.
pub struct MemoryTracker {
    current: AtomicUsize,
    limit: usize,
    peak: AtomicUsize,
}

impl MemoryTracker {
    #[must_use]
    pub fn new(limit: usize) -> Self {
        Self { current: AtomicUsize::new(0), limit, peak: AtomicUsize::new(0) }
    }

    /// Attempt to add `bytes` to the tracker without exceeding the limit.
    ///
    /// Returns `true` if the bytes were reserved, `false` if the reservation
    /// would exceed the global limit. This is the gate callers must use to
    /// enforce the global memory budget; `add()` is advisory-only and always
    /// succeeds.
    ///
    /// Uses a CAS loop: load current, reject if `current + bytes > limit`,
    /// otherwise `compare_exchange_weak`. Retries on CAS failure so concurrent
    /// producers race correctly — whichever producer wins the CAS reserves
    /// the bytes first.
    ///
    /// Peak tracking is updated on successful reservations.
    pub fn try_add(&self, bytes: usize) -> bool {
        let mut current = self.current.load(Ordering::Relaxed);
        let new = loop {
            let candidate = current.saturating_add(bytes);
            if candidate > self.limit {
                return false;
            }
            match self.current.compare_exchange_weak(
                current,
                candidate,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break candidate,
                Err(observed) => current = observed,
            }
        };
        // Update peak via CAS loop.
        let mut peak = self.peak.load(Ordering::Relaxed);
        while new > peak {
            match self.peak.compare_exchange_weak(peak, new, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => break,
                Err(observed) => peak = observed,
            }
        }
        true
    }

    /// Add `bytes` to the tracker.
    ///
    /// Always succeeds — backpressure is advisory. Callers check
    /// [`MemoryTracker::over_limit`] before producing work.
    ///
    /// Saturates at `usize::MAX` rather than overflowing, so callers that
    /// inadvertently pass absurd estimates will not corrupt the counter.
    /// Note that saturation still leaves the tracker in an "always over
    /// limit" state, which is the desired behaviour — it signals the
    /// pipeline to stop producing work.
    pub fn add(&self, bytes: usize) {
        // CAS loop for saturating add so the observed post-add value is
        // accurate even under concurrency.
        let mut current = self.current.load(Ordering::Relaxed);
        let new = loop {
            let candidate = current.saturating_add(bytes);
            match self.current.compare_exchange_weak(
                current,
                candidate,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break candidate,
                Err(observed) => current = observed,
            }
        };
        // Update peak via CAS loop.
        let mut peak = self.peak.load(Ordering::Relaxed);
        while new > peak {
            match self.peak.compare_exchange_weak(peak, new, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => break,
                Err(observed) => peak = observed,
            }
        }
    }

    /// Remove `bytes` from the tracker.
    ///
    /// Saturates at zero rather than wrapping, so a spurious double-sub
    /// will not silently produce a gigantic `current()` value.
    pub fn sub(&self, bytes: usize) {
        let mut current = self.current.load(Ordering::Relaxed);
        loop {
            let candidate = current.saturating_sub(bytes);
            match self.current.compare_exchange_weak(
                current,
                candidate,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(observed) => current = observed,
            }
        }
    }

    /// Current bytes tracked.
    #[must_use]
    pub fn current(&self) -> usize {
        self.current.load(Ordering::Relaxed)
    }

    /// Whether current bytes meet or exceed the limit.
    #[must_use]
    pub fn over_limit(&self) -> bool {
        self.current() >= self.limit
    }

    /// Whether current bytes are below half the limit (drain-mode exit hysteresis).
    #[must_use]
    pub fn drained(&self) -> bool {
        self.current() < self.limit / 2
    }

    /// Peak bytes observed.
    #[must_use]
    pub fn peak(&self) -> usize {
        self.peak.load(Ordering::Relaxed)
    }

    /// Limit in bytes.
    #[must_use]
    pub fn limit(&self) -> usize {
        self.limit
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_and_sub() {
        let mt = MemoryTracker::new(1000);
        mt.add(300);
        assert_eq!(mt.current(), 300);
        mt.add(200);
        assert_eq!(mt.current(), 500);
        mt.sub(100);
        assert_eq!(mt.current(), 400);
    }

    #[test]
    fn test_peak_tracking() {
        let mt = MemoryTracker::new(1000);
        mt.add(200);
        mt.add(300); // now 500
        mt.sub(100); // now 400
        mt.add(50); // now 450
        assert_eq!(mt.peak(), 500);
    }

    #[test]
    fn test_over_limit() {
        let mt = MemoryTracker::new(1000);
        assert!(!mt.over_limit());
        mt.add(999);
        assert!(!mt.over_limit());
        mt.add(1);
        assert!(mt.over_limit());
    }

    #[test]
    fn test_drained_hysteresis() {
        let mt = MemoryTracker::new(1000);
        mt.add(600);
        assert!(!mt.drained()); // at 600, limit/2 = 500, not drained
        mt.sub(200);
        assert!(mt.drained()); // at 400, below 500, drained
    }

    #[test]
    fn test_concurrent_adds() {
        use std::sync::Arc;
        use std::thread;

        let mt = Arc::new(MemoryTracker::new(10_000));
        let mut handles = vec![];
        for _ in 0..10 {
            let mt = mt.clone();
            handles.push(thread::spawn(move || {
                for _ in 0..100 {
                    mt.add(10);
                }
            }));
        }
        for h in handles {
            h.join().unwrap();
        }
        assert_eq!(mt.current(), 10_000);
        assert_eq!(mt.peak(), 10_000);
    }

    #[test]
    fn test_add_saturates_at_usize_max() {
        let mt = MemoryTracker::new(usize::MAX);
        mt.add(usize::MAX - 100);
        assert_eq!(mt.current(), usize::MAX - 100);
        // Adding more than remaining would overflow — should saturate.
        mt.add(1000);
        assert_eq!(mt.current(), usize::MAX);
        // Further adds stay at MAX.
        mt.add(1);
        assert_eq!(mt.current(), usize::MAX);
    }

    #[test]
    fn test_sub_saturates_at_zero() {
        let mt = MemoryTracker::new(1000);
        mt.add(50);
        // Spurious sub of more than is tracked would underflow — should saturate.
        mt.sub(1000);
        assert_eq!(mt.current(), 0);
    }

    #[test]
    fn test_try_add_succeeds_below_limit() {
        let mt = MemoryTracker::new(1000);
        assert!(mt.try_add(400));
        assert_eq!(mt.current(), 400);
        assert!(mt.try_add(500));
        assert_eq!(mt.current(), 900);
    }

    #[test]
    fn test_try_add_fails_at_limit() {
        let mt = MemoryTracker::new(1000);
        assert!(mt.try_add(900));
        // Would push us to 1100 > 1000 — must fail without mutating current.
        assert!(!mt.try_add(200));
        assert_eq!(mt.current(), 900);
        // Exactly at the limit is allowed.
        assert!(mt.try_add(100));
        assert_eq!(mt.current(), 1000);
        // Anything above zero now fails.
        assert!(!mt.try_add(1));
        assert_eq!(mt.current(), 1000);
    }

    #[test]
    fn test_try_add_updates_peak() {
        let mt = MemoryTracker::new(1000);
        assert!(mt.try_add(300));
        assert!(mt.try_add(200));
        assert_eq!(mt.peak(), 500);
        mt.sub(400);
        assert_eq!(mt.peak(), 500);
        assert!(mt.try_add(100));
        assert_eq!(mt.peak(), 500);
    }

    #[test]
    fn test_try_add_concurrent_reservations_are_bounded() {
        use std::sync::Arc;
        use std::sync::atomic::{AtomicUsize, Ordering as AO};
        use std::thread;

        // 10 threads each try 100 reservations of 10 bytes against a 500-byte
        // limit. The total reserved bytes must never exceed the limit; the
        // total number of successful reservations must be exactly limit/bytes.
        let mt = Arc::new(MemoryTracker::new(500));
        let successes = Arc::new(AtomicUsize::new(0));
        let mut handles = vec![];
        for _ in 0..10 {
            let mt = mt.clone();
            let successes = successes.clone();
            handles.push(thread::spawn(move || {
                for _ in 0..100 {
                    if mt.try_add(10) {
                        successes.fetch_add(1, AO::Relaxed);
                    }
                }
            }));
        }
        for h in handles {
            h.join().unwrap();
        }
        assert_eq!(mt.current(), 500);
        assert_eq!(successes.load(AO::Relaxed), 50);
    }
}

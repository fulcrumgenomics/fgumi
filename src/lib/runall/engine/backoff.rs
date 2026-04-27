//! Exponential backoff with jitter for hot spin-retry loops.
//!
//! A naive `std::thread::yield_now()` in a retry loop burns a core when the
//! producer/consumer cannot make progress. This module provides a small stateful
//! backoff that:
//! - Yields (no sleep syscall) on the first few attempts — cheap when the
//!   contention window is short.
//! - Grows exponentially from 10 µs to 1 ms on sustained contention — amortises
//!   the cost of repeated CAS failures or full-queue bounces.
//! - Applies ±25% jitter to every sleep — breaks up thundering-herd patterns
//!   where multiple workers wake up synchronously when an upstream queue drains.
//!
//! Callers use the pattern:
//! ```ignore
//! let mut backoff = Backoff::new();
//! loop {
//!     match try_thing() {
//!         Ok(v) => { backoff.reset(); return v; }
//!         Err(_) => backoff.snooze(),
//!     }
//! }
//! ```
//!
//! The constants match the v1 `unified_pipeline` tuning
//! (`MIN_BACKOFF_US=10`, `MAX_BACKOFF_US=1000`) which was settled after
//! iterated benchmarking on real workloads.

use std::time::Duration;

/// Number of `yield_now()` attempts before switching to sleeping.
///
/// Keep this small (default 4): the first few retries are likely to succeed
/// without a full syscall, but extending the yield-only phase too far wastes
/// CPU when the blocker is long-lived.
const YIELD_ATTEMPTS: u32 = 4;

/// Minimum sleep duration (microseconds).
const MIN_BACKOFF_US: u64 = 10;

/// Maximum sleep duration (microseconds).
const MAX_BACKOFF_US: u64 = 1_000;

/// Jitter factor (fraction of the current delay, applied as ±).
///
/// 0.25 = up to ±25% of the nominal delay is added/subtracted. Breaks up
/// synchronised wake-ups when many workers unblock at once.
const JITTER_FACTOR: f64 = 0.25;

/// Exponential backoff with jitter.
///
/// Stateful per-caller; not `Sync`. Each worker / retry-loop owns its own
/// instance so concurrent backoffs are independent.
#[derive(Debug)]
pub struct Backoff {
    /// Number of `snooze()` calls since the last `reset()`.
    attempts: u32,
    /// Current backoff delay in microseconds. Doubles on each sleep.
    current_us: u64,
}

impl Backoff {
    /// Create a new backoff with state reset.
    #[must_use]
    pub fn new() -> Self {
        Self { attempts: 0, current_us: MIN_BACKOFF_US }
    }

    /// Call when a retry is needed. The first few calls yield; subsequent
    /// calls sleep with exponentially increasing delay and ±25% jitter.
    pub fn snooze(&mut self) {
        if self.attempts < YIELD_ATTEMPTS {
            std::thread::yield_now();
            self.attempts += 1;
            return;
        }

        let jittered_us = Self::jitter(self.current_us);
        std::thread::sleep(Duration::from_micros(jittered_us));

        // Grow toward the cap. Using saturating_mul avoids overflow if the
        // cap is ever raised extreme values.
        self.current_us = self.current_us.saturating_mul(2).min(MAX_BACKOFF_US);
        self.attempts += 1;
    }

    /// Reset the backoff state. Call after a successful operation so the next
    /// stall starts at the minimum delay again.
    pub fn reset(&mut self) {
        self.attempts = 0;
        self.current_us = MIN_BACKOFF_US;
    }

    /// Apply ±`JITTER_FACTOR` random jitter to `us`. Pure for testability:
    /// the jittered value is always within `[us * (1 - JITTER_FACTOR),
    /// us * (1 + JITTER_FACTOR)]`.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::cast_precision_loss)]
    fn jitter(us: u64) -> u64 {
        use rand::RngExt;
        let delta = (us as f64) * JITTER_FACTOR;
        let offset = rand::rng().random_range(-delta..=delta);
        let jittered = (us as f64 + offset).max(0.0);
        jittered as u64
    }
}

impl Default for Backoff {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// First `YIELD_ATTEMPTS` snoozes should not grow `current_us` (they yield
    /// without sleeping). After that, each snooze doubles `current_us` up to
    /// the cap.
    #[test]
    fn test_backoff_monotonic_growth() {
        let mut b = Backoff::new();
        assert_eq!(b.current_us, MIN_BACKOFF_US);

        // Yield phase: current_us should not change.
        for _ in 0..YIELD_ATTEMPTS {
            b.snooze();
            assert_eq!(b.current_us, MIN_BACKOFF_US);
        }

        // Sleep phase: current_us doubles each call (capped).
        let mut last = b.current_us;
        for _ in 0..20 {
            b.snooze();
            assert!(b.current_us >= last, "backoff must be non-decreasing");
            assert!(b.current_us <= MAX_BACKOFF_US, "backoff must not exceed cap");
            last = b.current_us;
        }
        assert_eq!(b.current_us, MAX_BACKOFF_US);
    }

    #[test]
    fn test_backoff_reset_restores_initial_state() {
        let mut b = Backoff::new();
        for _ in 0..(YIELD_ATTEMPTS + 10) {
            b.snooze();
        }
        assert!(b.attempts > YIELD_ATTEMPTS);
        assert!(b.current_us > MIN_BACKOFF_US);

        b.reset();
        assert_eq!(b.attempts, 0);
        assert_eq!(b.current_us, MIN_BACKOFF_US);
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::cast_precision_loss)]
    fn test_jitter_bounds() {
        // Sample jitter many times and verify every result is within the
        // ±JITTER_FACTOR envelope.
        let base: u64 = 100;
        let delta = (base as f64) * JITTER_FACTOR;
        let lo = ((base as f64) - delta).max(0.0) as u64;
        let hi = ((base as f64) + delta) as u64;

        for _ in 0..1000 {
            let j = Backoff::jitter(base);
            assert!(j >= lo, "jittered {j} below lower bound {lo}");
            assert!(j <= hi, "jittered {j} above upper bound {hi}");
        }
    }

    #[test]
    fn test_default_matches_new() {
        let a = Backoff::new();
        let b = Backoff::default();
        assert_eq!(a.attempts, b.attempts);
        assert_eq!(a.current_us, b.current_us);
    }
}

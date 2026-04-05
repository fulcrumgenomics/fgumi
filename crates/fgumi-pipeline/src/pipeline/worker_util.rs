//! Shared utilities for worker threads: adaptive backoff and panic message extraction.

use std::time::Duration;

// ============================================================================
// WorkerBackoff
// ============================================================================

/// Adaptive backoff with per-worker jitter for work-stealing loops.
///
/// Progression: yield → 1µs → 2µs → 4µs → … → 100µs max.
/// Jitter is ±25% of the current backoff to reduce thundering-herd effects.
pub struct WorkerBackoff {
    backoff_us: u64,
    rng_state: u64,
}

impl WorkerBackoff {
    /// Create a new backoff seeded from the worker ID.
    #[must_use]
    pub fn new(worker_id: usize) -> Self {
        Self {
            backoff_us: 0,
            rng_state: (worker_id as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1,
        }
    }

    /// Reset backoff to zero on successful work.
    pub fn reset(&mut self) {
        self.backoff_us = 0;
    }

    /// Sleep for the current backoff duration with per-worker jitter.
    pub fn sleep(&mut self) {
        if self.backoff_us == 0 {
            std::thread::yield_now();
        } else {
            let jitter_range = self.backoff_us / 4;
            let jitter =
                if jitter_range > 0 { self.next_rng() % (jitter_range * 2 + 1) } else { 0 };
            let sleep_us = self.backoff_us.saturating_sub(jitter_range) + jitter;
            std::thread::sleep(Duration::from_micros(sleep_us));
        }
    }

    /// Increase backoff: yield → 1µs → 2µs → 4µs → … → 100µs max.
    pub fn increase(&mut self) {
        if self.backoff_us == 0 {
            self.backoff_us = 1;
        } else {
            self.backoff_us = (self.backoff_us * 2).min(100);
        }
    }

    /// Simple xorshift64 RNG.
    fn next_rng(&mut self) -> u64 {
        let mut x = self.rng_state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.rng_state = x;
        x
    }
}

// ============================================================================
// panic_message
// ============================================================================

/// Extract a human-readable message from a panic payload.
pub fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown panic".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_backoff_progression() {
        let mut b = WorkerBackoff::new(0);
        assert_eq!(b.backoff_us, 0);

        b.increase();
        assert_eq!(b.backoff_us, 1);

        b.increase();
        assert_eq!(b.backoff_us, 2);

        b.increase();
        assert_eq!(b.backoff_us, 4);

        for _ in 0..20 {
            b.increase();
        }
        assert_eq!(b.backoff_us, 100);

        b.reset();
        assert_eq!(b.backoff_us, 0);
    }

    #[test]
    fn test_panic_message() {
        let str_payload: Box<dyn std::any::Any + Send> = Box::new("test panic");
        assert_eq!(panic_message(&str_payload), "test panic");

        let string_payload: Box<dyn std::any::Any + Send> = Box::new(String::from("string panic"));
        assert_eq!(panic_message(&string_payload), "string panic");

        let other_payload: Box<dyn std::any::Any + Send> = Box::new(42i32);
        assert_eq!(panic_message(&other_payload), "unknown panic");
    }
}

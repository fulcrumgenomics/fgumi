//! `StepDrainCounter`: coordinates last-worker-wins for closing a step's
//! shared output queue when it reports `StepOutcome::Finished`.
//!
//! For `Parallel` steps, init to N (= worker count). Each clone returns
//! `Finished` independently when the shared input edge drains; each calls
//! `observe_drain`, and the clone that takes the counter to 0 is the "last
//! worker" — only it calls `mark_outputs_drained` (closing the shared output).
//! Otherwise a clone could close the output while a sibling is still pushing.
//!
//! For `Serial` / `Exclusive` steps, init to 1. The single finisher wins on
//! its first call.

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

#[derive(Debug)]
pub struct StepDrainCounter {
    remaining: AtomicUsize,
}

impl StepDrainCounter {
    /// Construct a counter with the given initial decrement budget. Returns
    /// an `Arc` for sharing across worker threads.
    #[must_use]
    pub fn new(initial: usize) -> Arc<Self> {
        Arc::new(Self { remaining: AtomicUsize::new(initial) })
    }

    /// Called by a worker when it observes drain on this step. Returns
    /// `true` if this is the last worker (counter went from 1 to 0).
    /// Subsequent calls (counter already 0) return `false`.
    pub fn observe_drain(&self) -> bool {
        // CAS-decrement loop: atomically decrement only if `> 0`. Avoids
        // underflow under concurrent over-calls.
        let mut prev = self.remaining.load(AtomicOrdering::Acquire);
        loop {
            if prev == 0 {
                return false;
            }
            match self.remaining.compare_exchange_weak(
                prev,
                prev - 1,
                AtomicOrdering::AcqRel,
                AtomicOrdering::Acquire,
            ) {
                Ok(_) => return prev == 1,
                Err(actual) => prev = actual,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_worker_with_init_one_wins() {
        let counter = StepDrainCounter::new(1);
        assert!(counter.observe_drain());
        assert!(!counter.observe_drain());
    }

    #[test]
    fn last_worker_with_init_n_wins() {
        let counter = StepDrainCounter::new(4);
        assert!(!counter.observe_drain());
        assert!(!counter.observe_drain());
        assert!(!counter.observe_drain());
        assert!(counter.observe_drain());
    }

    #[test]
    fn extra_calls_after_zero_return_false() {
        let counter = StepDrainCounter::new(2);
        counter.observe_drain();
        counter.observe_drain();
        assert!(!counter.observe_drain());
    }

    #[test]
    fn concurrent_decrements_have_exactly_one_winner() {
        use std::thread;
        let counter = StepDrainCounter::new(8);
        let winners: Vec<_> = (0..8)
            .map(|_| {
                let c = Arc::clone(&counter);
                thread::spawn(move || c.observe_drain())
            })
            .map(|h| h.join().unwrap())
            .collect();
        assert_eq!(winners.iter().filter(|&&w| w).count(), 1, "exactly one winner");
    }
}

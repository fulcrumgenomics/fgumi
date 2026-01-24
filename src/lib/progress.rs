//! Progress tracking utilities
//!
//! This module provides a thread-safe progress tracker for logging progress at regular intervals.
//! The tracker maintains an internal count and logs when interval boundaries are crossed.

use log::info;
use std::sync::atomic::{AtomicU64, Ordering};

/// Thread-safe progress tracker for logging progress at regular intervals.
///
/// Maintains an internal count and logs progress messages when the count crosses
/// interval boundaries. Safe to use from multiple threads.
///
/// # Example
/// ```
/// use fgumi_lib::progress::ProgressTracker;
///
/// let tracker = ProgressTracker::new("Processed records")
///     .with_interval(100);
///
/// // Add items and log at interval boundaries
/// for _ in 0..250 {
///     tracker.log_if_needed(1);  // Logs at 100, 200
/// }
/// tracker.log_final();  // Logs "Processed records 250 (complete)"
/// ```
///
/// # Multi-threaded Example
/// ```
/// use fgumi_lib::progress::ProgressTracker;
/// use std::sync::Arc;
///
/// let tracker = Arc::new(ProgressTracker::new("Processed records").with_interval(1000));
///
/// // Multiple threads can safely add to the same tracker
/// let tracker_clone = Arc::clone(&tracker);
/// std::thread::spawn(move || {
///     tracker_clone.log_if_needed(500);
/// });
/// ```
pub struct ProgressTracker {
    /// The logging interval - progress is logged when count crosses multiples of this.
    interval: u64,
    /// Message prefix for log output.
    message: String,
    /// Internal count of items processed (thread-safe).
    count: AtomicU64,
}

impl ProgressTracker {
    /// Create a new progress tracker with the specified message.
    ///
    /// The tracker starts with a count of 0 and a default interval of 10,000.
    ///
    /// # Arguments
    /// * `message` - Message prefix for progress logs (e.g., "Processed records")
    #[must_use]
    pub fn new(message: impl Into<String>) -> Self {
        Self { interval: 10_000, message: message.into(), count: AtomicU64::new(0) }
    }

    /// Set the logging interval.
    ///
    /// Progress will be logged each time the count crosses a multiple of this interval.
    /// For example, with interval=1000, logs will occur at 1000, 2000, 3000, etc.
    ///
    /// # Arguments
    /// * `interval` - The interval between progress logs
    #[must_use]
    pub fn with_interval(mut self, interval: u64) -> Self {
        self.interval = interval;
        self
    }

    /// Add to the count and log if an interval boundary was crossed.
    ///
    /// This method is thread-safe and can be called from multiple threads.
    /// It atomically adds `additional` to the internal count and logs progress
    /// for each interval boundary crossed.
    ///
    /// The behavior is equivalent to incrementing the counter one-by-one and
    /// checking at each step, but implemented efficiently with a single atomic
    /// operation.
    ///
    /// # Arguments
    /// * `additional` - Number of items to add to the count
    ///
    /// # Returns
    /// `true` if the final count is exactly a multiple of the interval,
    /// `false` otherwise. This is useful for `log_final()` to know if a
    /// final message is needed.
    ///
    /// # Example
    /// ```
    /// use fgumi_lib::progress::ProgressTracker;
    ///
    /// let tracker = ProgressTracker::new("Items").with_interval(100);
    /// tracker.log_if_needed(50);   // count=50, returns false, no log
    /// tracker.log_if_needed(60);   // count=110, returns false, logs "Items 100"
    /// tracker.log_if_needed(90);   // count=200, returns true, logs "Items 200"
    /// ```
    pub fn log_if_needed(&self, additional: u64) -> bool {
        if additional == 0 {
            // No change, just check if current count is on interval
            let count = self.count.load(Ordering::Relaxed);
            return count > 0 && count.is_multiple_of(self.interval);
        }

        let prev = self.count.fetch_add(additional, Ordering::Relaxed);
        let new_count = prev + additional;

        // Calculate how many interval boundaries we crossed
        let prev_intervals = prev / self.interval;
        let new_intervals = new_count / self.interval;

        // Log for each interval boundary crossed
        for i in (prev_intervals + 1)..=new_intervals {
            let milestone = i * self.interval;
            info!("{} {}", self.message, milestone);
        }

        // Return true if we landed exactly on an interval
        new_count.is_multiple_of(self.interval)
    }

    /// Log final progress.
    ///
    /// If the current count is not exactly on an interval boundary (i.e., if
    /// `log_if_needed(0)` returns `false`), logs a final message with "(complete)".
    /// If the count is exactly on an interval, the last `log_if_needed` call
    /// already logged it, so no additional message is needed.
    ///
    /// # Example
    /// ```
    /// use fgumi_lib::progress::ProgressTracker;
    ///
    /// let tracker = ProgressTracker::new("Items").with_interval(100);
    /// tracker.log_if_needed(250);  // Logs "Items 100", "Items 200"
    /// tracker.log_final();          // Logs "Items 250 (complete)"
    ///
    /// let tracker2 = ProgressTracker::new("Items").with_interval(100);
    /// tracker2.log_if_needed(200);  // Logs "Items 100", "Items 200"
    /// tracker2.log_final();          // No log (200 is exactly on interval)
    /// ```
    pub fn log_final(&self) {
        if !self.log_if_needed(0) {
            let count = self.count.load(Ordering::Relaxed);
            if count > 0 {
                info!("{} {} (complete)", self.message, count);
            }
        }
    }

    /// Get the current count.
    ///
    /// # Returns
    /// The current count of items processed.
    #[must_use]
    pub fn count(&self) -> u64 {
        self.count.load(Ordering::Relaxed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_progress_tracker_new() {
        let tracker = ProgressTracker::new("Processing");
        assert_eq!(tracker.interval, 10_000);
        assert_eq!(tracker.message, "Processing");
        assert_eq!(tracker.count(), 0);
    }

    #[test]
    fn test_progress_tracker_with_interval() {
        let tracker = ProgressTracker::new("Processing").with_interval(100);
        assert_eq!(tracker.interval, 100);
    }

    #[test]
    fn test_log_if_needed_returns_correctly() {
        let tracker = ProgressTracker::new("Test").with_interval(10);

        // Not on interval
        assert!(!tracker.log_if_needed(5)); // count=5
        assert!(!tracker.log_if_needed(3)); // count=8

        // Crosses interval, lands on it
        assert!(tracker.log_if_needed(2)); // count=10, exactly on interval

        // Not on interval
        assert!(!tracker.log_if_needed(5)); // count=15

        // Crosses interval, doesn't land on it
        assert!(!tracker.log_if_needed(10)); // count=25, crossed 20
    }

    #[test]
    fn test_log_if_needed_zero() {
        let tracker = ProgressTracker::new("Test").with_interval(10);

        // Zero count, zero additional
        assert!(!tracker.log_if_needed(0));

        // Add to exactly on interval
        tracker.log_if_needed(10);
        assert!(tracker.log_if_needed(0)); // count=10, exactly on interval

        // Add more, not on interval
        tracker.log_if_needed(5);
        assert!(!tracker.log_if_needed(0)); // count=15
    }

    #[test]
    fn test_count() {
        let tracker = ProgressTracker::new("Test").with_interval(100);

        assert_eq!(tracker.count(), 0);
        tracker.log_if_needed(50);
        assert_eq!(tracker.count(), 50);
        tracker.log_if_needed(75);
        assert_eq!(tracker.count(), 125);
    }

    #[test]
    fn test_crossing_multiple_intervals() {
        let tracker = ProgressTracker::new("Test").with_interval(10);

        // Cross multiple intervals at once (10, 20, 30)
        assert!(!tracker.log_if_needed(35)); // count=35, crossed 10, 20, 30 but not on interval
        assert_eq!(tracker.count(), 35);

        // Cross to exactly on interval
        assert!(tracker.log_if_needed(5)); // count=40
    }

    #[test]
    fn test_thread_safety() {
        use std::sync::Arc;
        use std::thread;

        let tracker = Arc::new(ProgressTracker::new("Test").with_interval(1000));
        let mut handles = vec![];

        // Spawn 10 threads, each adding 100 items
        for _ in 0..10 {
            let tracker_clone = Arc::clone(&tracker);
            handles.push(thread::spawn(move || {
                for _ in 0..100 {
                    tracker_clone.log_if_needed(1);
                }
            }));
        }

        for handle in handles {
            handle.join().unwrap();
        }

        // Total should be 1000
        assert_eq!(tracker.count(), 1000);
    }
}

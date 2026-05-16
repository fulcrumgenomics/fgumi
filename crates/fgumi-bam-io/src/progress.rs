//! Progress tracking utilities
//!
//! This module provides a thread-safe progress tracker for logging progress at regular intervals.
//! The tracker maintains an internal count and logs when interval boundaries are crossed.
//! When a total is known, it also displays percentage complete and ETA using an exponential
//! moving average (EMA) of the processing rate with bias correction (tqdm-style).

use log::info;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;
use std::time::Instant;

/// Smoothing constant for the EMA rate estimator.
/// 0.3 balances responsiveness to rate changes with stability.
/// Same default as tqdm.
const EMA_ALPHA: f64 = 0.3;

/// State for the exponential moving average rate estimator.
struct EmaState {
    /// Smoothed rate (records per second), pre-bias-correction.
    smoothed_rate: f64,
    /// Number of EMA updates (for bias correction).
    calls: u32,
    /// Count at last EMA update.
    last_count: u64,
    /// Time at last EMA update.
    last_time: Instant,
}

impl EmaState {
    fn new() -> Self {
        Self { smoothed_rate: 0.0, calls: 0, last_count: 0, last_time: Instant::now() }
    }

    /// Update the EMA with a new observation and return the bias-corrected rate.
    fn update(&mut self, current_count: u64) -> f64 {
        if current_count <= self.last_count {
            return self.corrected_rate();
        }

        let now = Instant::now();
        let dt = now.duration_since(self.last_time).as_secs_f64();
        if dt > 0.0 {
            #[allow(clippy::cast_precision_loss)]
            let dn = (current_count - self.last_count) as f64;
            let instantaneous_rate = dn / dt;
            self.smoothed_rate =
                EMA_ALPHA * instantaneous_rate + (1.0 - EMA_ALPHA) * self.smoothed_rate;
            self.calls += 1;
            self.last_count = current_count;
            self.last_time = now;
        }
        self.corrected_rate()
    }

    /// Return the bias-corrected rate estimate.
    ///
    /// Uses the correction factor `1 / (1 - (1-α)^n)` to compensate for
    /// zero-initialization of the EMA, giving accurate estimates even with
    /// only a few updates.
    fn corrected_rate(&self) -> f64 {
        if self.calls == 0 {
            return 0.0;
        }
        let beta = 1.0 - EMA_ALPHA;
        let correction = 1.0 - beta.powi(self.calls.cast_signed());
        if correction <= 0.0 { 0.0 } else { self.smoothed_rate / correction }
    }
}

/// Formats a [`Duration`] as a human-readable string (e.g. `"2m 15s"`, `"1h 30m"`).
fn format_duration(duration: Duration) -> String {
    let secs = duration.as_secs();
    if secs < 60 {
        format!("{secs}s")
    } else if secs < 3600 {
        let mins = secs / 60;
        let remaining_secs = secs % 60;
        if remaining_secs == 0 { format!("{mins}m") } else { format!("{mins}m {remaining_secs}s") }
    } else {
        let hours = secs / 3600;
        let mins = (secs % 3600) / 60;
        if mins == 0 { format!("{hours}h") } else { format!("{hours}h {mins}m") }
    }
}

/// Convert seconds (f64) to a formatted duration string.
fn fmt_duration(secs: f64) -> String {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    format_duration(Duration::from_secs(secs.round() as u64))
}

/// Thread-safe progress tracker for logging progress at regular intervals.
///
/// Maintains an internal count and logs progress messages when the count crosses
/// interval boundaries. Safe to use from multiple threads.
///
/// When a total is set via [`with_total`](Self::with_total), progress messages include
/// percentage complete and an ETA estimate using an exponential moving average of the
/// processing rate with bias correction.
///
/// # Example
/// ```
/// use fgumi_bam_io::progress::ProgressTracker;
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
/// use fgumi_bam_io::progress::ProgressTracker;
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
    /// Next count value at which a milestone should be logged. Monotonically
    /// advanced by the slow path under the EMA mutex. Reading this on the fast
    /// path (instead of dividing `count` by `interval`) cuts ~30-60 cycles per
    /// call — at ~32 M `log_if_needed(1)` calls/s on a 16-thread sort that
    /// shows up as ~1 % of total CPU in profiles.
    next_threshold: AtomicU64,
    /// Optional total count for percentage and ETA display.
    total: Option<u64>,
    /// Time the tracker was created (for elapsed time in final message).
    start_time: Instant,
    /// EMA rate estimator state (only accessed during logging, so contention is negligible).
    ema: Mutex<EmaState>,
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
        let interval = 10_000;
        Self {
            interval,
            message: message.into(),
            count: AtomicU64::new(0),
            next_threshold: AtomicU64::new(interval),
            total: None,
            start_time: Instant::now(),
            ema: Mutex::new(EmaState::new()),
        }
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
        // `with_interval` is part of the builder chain (single-threaded), so a
        // non-atomic write to `next_threshold` would be sound here, but a
        // `store` is just as cheap and keeps the field's invariant local.
        self.next_threshold.store(interval, Ordering::Relaxed);
        self
    }

    /// Set the total expected count.
    ///
    /// When set, progress messages include percentage complete and an ETA estimate
    /// using an exponential moving average of the processing rate.
    ///
    /// # Arguments
    /// * `total` - The total expected count of items
    #[must_use]
    pub fn with_total(mut self, total: u64) -> Self {
        self.total = (total > 0).then_some(total);
        self
    }

    /// Add to the count and log if an interval boundary was crossed.
    ///
    /// This method is thread-safe and can be called from multiple threads.
    /// It atomically adds `additional` to the internal count and logs progress
    /// for each interval boundary crossed.
    ///
    /// When a total is set, log messages include percentage and ETA.
    ///
    /// # Arguments
    /// * `additional` - Number of items to add to the count
    ///
    /// # Returns
    /// `true` if the final count is exactly a multiple of the interval,
    /// `false` otherwise. This is useful for `log_final()` to know if a
    /// final message is needed.
    pub fn log_if_needed(&self, additional: u64) -> bool {
        if additional == 0 {
            // Cold path: only called from `log_final`. The modulo is fine here.
            let count = self.count.load(Ordering::Relaxed);
            return count > 0 && count.is_multiple_of(self.interval);
        }

        let prev = self.count.fetch_add(additional, Ordering::Relaxed);
        let new_count = prev + additional;

        // Fast path: count is monotonically non-decreasing and the interval
        // is fixed, so the thresholds form an arithmetic sequence (interval,
        // 2*interval, ...). We only need to know the next one and compare.
        // The common case is "haven't crossed yet" — predicted not-taken.
        let threshold = self.next_threshold.load(Ordering::Relaxed);
        if new_count < threshold {
            return false;
        }

        self.handle_crossing(new_count)
    }

    /// Slow path: log every milestone we just crossed, advance `next_threshold`.
    #[cold]
    #[allow(clippy::cast_precision_loss)]
    fn handle_crossing(&self, new_count: u64) -> bool {
        // Take the EMA lock first — it serializes both the EMA update and the
        // `next_threshold` advance, so concurrent crossings can't double-log.
        let Ok(mut ema) = self.ema.lock() else {
            // Mutex poisoned — give up on logging this call. Should never
            // happen in practice since the slow path doesn't panic.
            return false;
        };

        // Re-read under the lock: another thread may have advanced past us
        // while we were waiting. If so, our milestones were already logged.
        let mut threshold = self.next_threshold.load(Ordering::Relaxed);
        if new_count < threshold {
            return false;
        }

        let rate = if self.total.is_some() { ema.update(new_count) } else { 0.0 };

        while threshold <= new_count {
            let milestone = threshold;
            if let Some(total) = self.total {
                let pct = (milestone as f64 / total as f64) * 100.0;
                // Derive remaining work from milestone, not new_count, so each
                // logged line shows the ETA appropriate for that milestone.
                let eta_suffix = if rate > 0.0 {
                    let remaining = total.saturating_sub(milestone) as f64;
                    format!(", ETA ~{}", fmt_duration(remaining / rate))
                } else {
                    String::new()
                };
                info!("{} {} / {} ({:.1}%{})", self.message, milestone, total, pct, eta_suffix);
            } else {
                info!("{} {}", self.message, milestone);
            }
            threshold += self.interval;
        }
        self.next_threshold.store(threshold, Ordering::Relaxed);

        // True iff `new_count` landed exactly on the most recently logged
        // milestone, which is `threshold - self.interval` after the advance.
        new_count + self.interval == threshold
    }

    /// Log final progress.
    ///
    /// When a total is set, always logs a completion message with elapsed time.
    /// Otherwise, logs only if the current count is not on an interval boundary.
    pub fn log_final(&self) {
        let count = self.count.load(Ordering::Relaxed);
        if count == 0 && self.total.is_none() {
            return;
        }

        if self.total.is_some() {
            let elapsed = self.start_time.elapsed().as_secs_f64();
            info!("{} {} (complete, {})", self.message, count, fmt_duration(elapsed));
        } else if !self.log_if_needed(0) {
            info!("{} {} (complete)", self.message, count);
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
    use rstest::rstest;

    use super::*;

    #[test]
    fn test_progress_tracker_new() {
        let tracker = ProgressTracker::new("Processing");
        assert_eq!(tracker.interval, 10_000);
        assert_eq!(tracker.message, "Processing");
        assert_eq!(tracker.count(), 0);
        assert!(tracker.total.is_none());
    }

    #[test]
    fn test_progress_tracker_with_interval() {
        let tracker = ProgressTracker::new("Processing").with_interval(100);
        assert_eq!(tracker.interval, 100);
    }

    #[test]
    fn test_progress_tracker_with_total() {
        let tracker = ProgressTracker::new("Processing").with_total(1000);
        assert_eq!(tracker.total, Some(1000));
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
            handle.join().expect("thread should join successfully");
        }

        // Total should be 1000
        assert_eq!(tracker.count(), 1000);
    }

    #[test]
    fn test_with_total_tracks_count() {
        let tracker = ProgressTracker::new("Test").with_interval(10).with_total(100);

        tracker.log_if_needed(25);
        assert_eq!(tracker.count(), 25);
        tracker.log_if_needed(75);
        assert_eq!(tracker.count(), 100);
    }

    #[rstest]
    #[case(0.0, "0s")]
    #[case(59.0, "59s")]
    #[case(59.5, "1m")]
    #[case(90.0, "1m 30s")]
    #[case(3600.0, "1h")]
    #[case(5400.0, "1h 30m")]
    fn test_fmt_duration(#[case] secs: f64, #[case] expected: &str) {
        assert_eq!(fmt_duration(secs), expected);
    }

    #[test]
    fn test_ema_bias_correction() {
        let mut ema = EmaState::new();

        // With zero calls, corrected rate should be 0
        assert!(ema.corrected_rate().abs() < f64::EPSILON);

        // After first update, corrected rate equals instantaneous rate
        // (bias correction factor is 1/(1-0.7^1) = 1/0.3 = 3.33,
        //  and smoothed_rate = 0.3 * rate, so corrected = rate)
        std::thread::sleep(std::time::Duration::from_millis(10));
        ema.last_count = 0;
        let rate = ema.update(1000);
        assert!(rate > 0.0, "rate should be positive after first update");
    }
}

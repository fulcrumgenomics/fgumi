//! Per-stage performance counters for the parallel pipeline.

use std::sync::atomic::{AtomicU64, Ordering};

/// Per-stage performance counters for the parallel pipeline.
///
/// Thread-safe via atomic operations. Use [`StageStats::record`] from any worker thread
/// to accumulate timing and item counts, then read totals after the pipeline completes.
pub struct StageStats {
    /// Total nanoseconds spent processing items in this stage.
    total_nanos: AtomicU64,
    /// Number of items processed.
    items_processed: AtomicU64,
    /// Stage name for logging.
    name: &'static str,
}

impl StageStats {
    /// Create a new stats tracker for the named stage.
    #[must_use]
    pub fn new(name: &'static str) -> Self {
        Self { total_nanos: AtomicU64::new(0), items_processed: AtomicU64::new(0), name }
    }

    /// Record the elapsed time and item count for one `process()` call.
    pub fn record(&self, nanos: u64, items: u64) {
        self.total_nanos.fetch_add(nanos, Ordering::Relaxed);
        self.items_processed.fetch_add(items, Ordering::Relaxed);
    }

    /// Total wall-clock seconds accumulated across all worker threads.
    #[must_use]
    #[expect(clippy::cast_precision_loss, reason = "nanosecond precision loss is acceptable")]
    pub fn total_secs(&self) -> f64 {
        self.total_nanos.load(Ordering::Relaxed) as f64 / 1_000_000_000.0
    }

    /// Total number of items processed.
    #[must_use]
    pub fn items(&self) -> u64 {
        self.items_processed.load(Ordering::Relaxed)
    }

    /// Stage name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        self.name
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[expect(clippy::float_cmp, reason = "exact float comparison is fine for test values")]
    fn test_record_and_read() {
        let stats = StageStats::new("test_stage");
        assert_eq!(stats.items(), 0);
        assert_eq!(stats.total_secs(), 0.0);
        assert_eq!(stats.name(), "test_stage");

        stats.record(1_000_000_000, 10);
        assert_eq!(stats.items(), 10);
        assert!((stats.total_secs() - 1.0).abs() < f64::EPSILON);

        stats.record(500_000_000, 5);
        assert_eq!(stats.items(), 15);
        assert!((stats.total_secs() - 1.5).abs() < f64::EPSILON);
    }
}

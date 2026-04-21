//! Backpressure sampling: summarizes queue and memory state for the scheduler.

use std::sync::Arc;

use super::memory::MemoryTracker;

const HIGH_WATERMARK: f64 = 0.75;
const LOW_WATERMARK: f64 = 0.25;

/// Snapshot of pipeline pressure, sampled once per scheduler iteration.
#[derive(Debug, Clone)]
pub struct BackpressureState {
    /// Per-queue fill summary. Index = queue ID (one per inter-stage queue).
    pub queue_summaries: Vec<QueueSummary>,
    /// Whether global memory is at or above limit.
    pub memory_high: bool,
    /// Whether global memory has drained to the hysteresis threshold.
    pub memory_drained: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct QueueSummary {
    /// Queue ID (index in the pipeline's queue vector).
    pub id: usize,
    /// Slot fill ratio (0.0 to 1.0).
    pub slot_fill: f64,
    /// Memory fill ratio (0.0 to 1.0+).
    pub memory_fill: f64,
    /// Whether the queue's producer has closed it.
    pub closed: bool,
    /// Whether the queue is empty.
    pub empty: bool,
}

impl QueueSummary {
    /// Whether this queue is at high-water mark (under pressure).
    #[must_use]
    pub fn is_high(&self) -> bool {
        self.slot_fill >= HIGH_WATERMARK || self.memory_fill >= HIGH_WATERMARK
    }

    /// Whether this queue is at low-water mark (backpressure can drain).
    #[must_use]
    pub fn is_low(&self) -> bool {
        self.slot_fill <= LOW_WATERMARK && self.memory_fill <= LOW_WATERMARK
    }
}

impl BackpressureState {
    /// Construct from per-queue summaries and global memory state.
    pub fn new(queue_summaries: Vec<QueueSummary>, global_memory: &Arc<MemoryTracker>) -> Self {
        Self {
            queue_summaries,
            memory_high: global_memory.over_limit(),
            memory_drained: global_memory.drained(),
        }
    }

    /// Is any queue under pressure (high-water mark)?
    #[must_use]
    pub fn any_high(&self) -> bool {
        self.queue_summaries.iter().any(QueueSummary::is_high)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_queue_summary_high_low() {
        let high =
            QueueSummary { id: 0, slot_fill: 0.9, memory_fill: 0.5, closed: false, empty: false };
        assert!(high.is_high());
        assert!(!high.is_low());

        let low =
            QueueSummary { id: 0, slot_fill: 0.1, memory_fill: 0.1, closed: false, empty: false };
        assert!(!low.is_high());
        assert!(low.is_low());

        let mid =
            QueueSummary { id: 0, slot_fill: 0.5, memory_fill: 0.5, closed: false, empty: false };
        assert!(!mid.is_high());
        assert!(!mid.is_low());
    }

    #[test]
    fn test_backpressure_memory_flags() {
        let tracker = Arc::new(MemoryTracker::new(100));
        tracker.add(50);
        let bp = BackpressureState::new(vec![], &tracker);
        assert!(!bp.memory_high);
        assert!(!bp.memory_drained); // 50 >= 50 (half of 100)

        tracker.add(60); // now 110, over limit
        let bp = BackpressureState::new(vec![], &tracker);
        assert!(bp.memory_high);
    }

    #[test]
    fn test_backpressure_any_high() {
        let tracker = Arc::new(MemoryTracker::new(1000));
        let low =
            QueueSummary { id: 0, slot_fill: 0.1, memory_fill: 0.1, closed: false, empty: false };
        let high =
            QueueSummary { id: 1, slot_fill: 0.9, memory_fill: 0.5, closed: false, empty: false };
        let bp = BackpressureState::new(vec![low, high], &tracker);
        assert!(bp.any_high());

        let bp2 = BackpressureState::new(vec![low], &tracker);
        assert!(!bp2.any_high());
    }
}

//! Watchdog: detects deadlocks by monitoring queue progress.
//!
//! The watchdog runs in its own thread. It samples queue push/pop counters
//! each interval and reports a deadlock if no counter has changed for the
//! configured timeout.
//!
//! Each [`QueueProgress`] also records the most recent push and pop
//! timestamps (milliseconds since the tracker's start instant) so that on a
//! deadlock the watchdog can attribute the stall to a specific queue:
//! "queue 3: last push 12s ago, last pop 3s ago" points directly at the
//! stuck stage.

use std::sync::Arc;
use std::sync::LazyLock;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use super::cancel::CancelToken;

/// Process-global start instant used as the reference for push/pop
/// timestamps. Lazily initialised on first access. Using a single global
/// instant (rather than one `Instant` per `QueueProgress`) keeps the
/// structure `Default`-constructible and lets timestamps from different
/// queues be compared directly.
static PROGRESS_EPOCH: LazyLock<Instant> = LazyLock::new(Instant::now);

/// Milliseconds since [`PROGRESS_EPOCH`]. Saturates to `u64::MAX` if the
/// process ever runs long enough to overflow (a few hundred million years).
///
/// Clamps the minimum returned value to 1: we reserve 0 as the "no events
/// recorded" sentinel in [`QueueProgress`], so callers can use
/// `load() == 0` to distinguish never-touched from just-touched.
#[allow(clippy::cast_possible_truncation)]
fn now_millis() -> u64 {
    let elapsed = PROGRESS_EPOCH.elapsed();
    let ms = elapsed.as_millis();
    let clamped = if ms > u128::from(u64::MAX) { u64::MAX } else { ms as u64 };
    clamped.max(1)
}

/// Counters tracked per queue for deadlock detection.
#[derive(Default)]
pub struct QueueProgress {
    pub pushes: AtomicU64,
    pub pops: AtomicU64,
    /// Milliseconds since [`PROGRESS_EPOCH`] of the most recent successful
    /// push. Zero means "no pushes yet"; callers use
    /// [`QueueProgress::snapshot_diagnostics`] to format human-readable
    /// ages.
    pub last_push_ms: AtomicU64,
    /// Milliseconds since [`PROGRESS_EPOCH`] of the most recent successful
    /// pop.
    pub last_pop_ms: AtomicU64,
}

impl QueueProgress {
    pub fn note_push(&self) {
        self.pushes.fetch_add(1, Ordering::Relaxed);
        self.last_push_ms.store(now_millis(), Ordering::Relaxed);
    }

    pub fn note_pop(&self) {
        self.pops.fetch_add(1, Ordering::Relaxed);
        self.last_pop_ms.store(now_millis(), Ordering::Relaxed);
    }

    #[must_use]
    fn snapshot(&self) -> (u64, u64) {
        (self.pushes.load(Ordering::Relaxed), self.pops.load(Ordering::Relaxed))
    }

    /// Format a short human-readable diagnostic line for this queue. Called
    /// by the watchdog on deadlock so the error message points to the stuck
    /// queue rather than a generic "pipeline deadlocked".
    #[must_use]
    pub fn snapshot_diagnostics(&self, queue_index: usize) -> String {
        let now = now_millis();
        let last_push = self.last_push_ms.load(Ordering::Relaxed);
        let last_pop = self.last_pop_ms.load(Ordering::Relaxed);
        let pushes = self.pushes.load(Ordering::Relaxed);
        let pops = self.pops.load(Ordering::Relaxed);
        let push_age = Self::format_age_ms(last_push, now);
        let pop_age = Self::format_age_ms(last_pop, now);
        format!(
            "queue {queue_index}: pushes={pushes} pops={pops} last_push={push_age} last_pop={pop_age}"
        )
    }

    /// Format an age as `never` (when the timestamp is zero, i.e. no events
    /// have been recorded yet) or `<N>ms ago`. Saturates on any unexpected
    /// clock regression so the caller never sees a negative value.
    fn format_age_ms(last_ms: u64, now_ms: u64) -> String {
        if last_ms == 0 {
            "never".to_string()
        } else {
            let age = now_ms.saturating_sub(last_ms);
            format!("{age}ms ago")
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct WatchdogConfig {
    pub timeout: Duration,
    pub poll_interval: Duration,
}

impl Default for WatchdogConfig {
    fn default() -> Self {
        Self { timeout: Duration::from_secs(10), poll_interval: Duration::from_millis(500) }
    }
}

/// Result of running the watchdog once.
#[derive(Debug, PartialEq, Eq)]
pub enum WatchdogOutcome {
    Ok,
    Deadlock,
    Cancelled,
}

/// Monitor queue progress until cancelled or deadlock detected.
///
/// Blocks the calling thread. Typically run in a dedicated watchdog thread.
/// On deadlock, emits a `tracing::warn!` line per queue with the most recent
/// push/pop ages so operators can attribute the stall to a specific stage.
#[must_use]
pub fn watch(
    progress: &[Arc<QueueProgress>],
    cancel: &CancelToken,
    cfg: WatchdogConfig,
) -> WatchdogOutcome {
    let mut last_snapshots: Vec<(u64, u64)> = progress.iter().map(|p| p.snapshot()).collect();
    let mut last_change = Instant::now();

    loop {
        if cancel.is_cancelled() {
            return WatchdogOutcome::Cancelled;
        }

        std::thread::sleep(cfg.poll_interval);

        let now_snapshots: Vec<(u64, u64)> = progress.iter().map(|p| p.snapshot()).collect();
        if now_snapshots != last_snapshots {
            last_change = Instant::now();
            last_snapshots = now_snapshots;
        } else if last_change.elapsed() >= cfg.timeout {
            tracing::warn!(
                "pipeline deadlock suspected: no queue progress for {}s",
                cfg.timeout.as_secs()
            );
            for (i, p) in progress.iter().enumerate() {
                tracing::warn!("pipeline deadlock diagnostic: {}", p.snapshot_diagnostics(i));
            }
            return WatchdogOutcome::Deadlock;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_watchdog_ok_when_progress_continues() {
        let progress = vec![Arc::new(QueueProgress::default())];
        let progress_for_producer = progress[0].clone();
        let cancel = CancelToken::new();

        let cancel_for_producer = cancel.clone();
        let producer = std::thread::spawn(move || {
            for _ in 0..10 {
                progress_for_producer.note_push();
                std::thread::sleep(Duration::from_millis(50));
            }
            cancel_for_producer.cancel();
        });

        let outcome = watch(
            &progress,
            &cancel,
            WatchdogConfig {
                timeout: Duration::from_secs(5),
                poll_interval: Duration::from_millis(20),
            },
        );

        producer.join().unwrap();
        assert_eq!(outcome, WatchdogOutcome::Cancelled);
    }

    #[test]
    fn test_watchdog_detects_deadlock() {
        let progress = vec![Arc::new(QueueProgress::default())];
        let cancel = CancelToken::new();

        let outcome = watch(
            &progress,
            &cancel,
            WatchdogConfig {
                timeout: Duration::from_millis(200),
                poll_interval: Duration::from_millis(20),
            },
        );
        assert_eq!(outcome, WatchdogOutcome::Deadlock);
    }

    #[test]
    fn test_watchdog_respects_cancellation() {
        let progress = vec![Arc::new(QueueProgress::default())];
        let cancel = CancelToken::new();
        cancel.cancel();

        let outcome = watch(
            &progress,
            &cancel,
            WatchdogConfig {
                timeout: Duration::from_secs(10),
                poll_interval: Duration::from_millis(20),
            },
        );
        assert_eq!(outcome, WatchdogOutcome::Cancelled);
    }

    #[test]
    fn test_push_pop_update_timestamps_monotonically() {
        let p = QueueProgress::default();
        assert_eq!(p.last_push_ms.load(Ordering::Relaxed), 0);
        assert_eq!(p.last_pop_ms.load(Ordering::Relaxed), 0);

        p.note_push();
        let first_push = p.last_push_ms.load(Ordering::Relaxed);
        assert!(first_push > 0, "first push should set a non-zero timestamp");

        // Small sleep so the next timestamp has a chance to be larger.
        std::thread::sleep(Duration::from_millis(5));

        p.note_push();
        let second_push = p.last_push_ms.load(Ordering::Relaxed);
        assert!(
            second_push >= first_push,
            "second push timestamp {second_push} should be >= first {first_push}"
        );

        p.note_pop();
        let first_pop = p.last_pop_ms.load(Ordering::Relaxed);
        assert!(first_pop > 0);
    }

    #[test]
    fn test_snapshot_diagnostics_initial_state() {
        let p = QueueProgress::default();
        let s = p.snapshot_diagnostics(3);
        assert!(s.contains("queue 3"));
        assert!(s.contains("pushes=0"));
        assert!(s.contains("pops=0"));
        // No events yet => "never" for both timestamps.
        assert!(s.contains("last_push=never"));
        assert!(s.contains("last_pop=never"));
    }

    #[test]
    fn test_snapshot_diagnostics_after_push() {
        let p = QueueProgress::default();
        p.note_push();
        p.note_push();
        let s = p.snapshot_diagnostics(0);
        assert!(s.contains("queue 0"));
        assert!(s.contains("pushes=2"));
        assert!(s.contains("pops=0"));
        assert!(s.contains("ms ago"));
        assert!(s.contains("last_pop=never"));
    }
}

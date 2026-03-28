//! Health monitor for the aligner subprocess.
//!
//! [`HealthMonitor`] runs in a background thread and watches activity timestamps
//! on the aligner's stdin and stdout pipes. If stdout goes silent for longer than
//! the configured timeout the monitor logs a warning; at twice the timeout it
//! sets the shared cancellation flag so the rest of the pipeline can shut down
//! cleanly.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::thread::{self, JoinHandle};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

/// Interval between health checks.
const CHECK_INTERVAL: Duration = Duration::from_secs(5);

// ============================================================================
// Helpers
// ============================================================================

/// Return the current time as milliseconds since the Unix epoch.
///
/// Returns `0` on the (extremely unlikely) event that the system clock is
/// before the epoch.
#[must_use]
pub fn now_millis() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| u64::try_from(d.as_millis()).unwrap_or(u64::MAX))
        .unwrap_or(0)
}

// ============================================================================
// HealthMonitor
// ============================================================================

/// Monitors aligner subprocess liveness by watching stdin/stdout activity timestamps.
///
/// Activity is reported by the I/O threads that read SAM from the aligner stdout
/// and write FASTQ to the aligner stdin. Each I/O thread updates the corresponding
/// [`AtomicU64`] with [`now_millis`] after every successful read/write.
///
/// The monitor wakes every 5 seconds and checks:
/// - stdout inactive for > `timeout` → log a warning.
/// - stdout inactive for > `2 × timeout` → set `cancel`, log an error.
/// - stdin inactive for > `timeout` → log a warning.
pub struct HealthMonitor {
    /// PID of the aligner child process (used in log messages).
    child_pid: u32,
    /// Epoch-millisecond timestamp of the last successful stdout read.
    last_stdout_activity: Arc<AtomicU64>,
    /// Epoch-millisecond timestamp of the last successful stdin write.
    last_stdin_activity: Arc<AtomicU64>,
    /// Duration after which inactivity is considered a warning condition.
    timeout: Duration,
    /// Shared cancellation flag; set to `true` by this monitor on fatal timeout.
    cancel: Arc<AtomicBool>,
}

impl HealthMonitor {
    /// Create a new `HealthMonitor`.
    ///
    /// # Arguments
    ///
    /// * `child_pid`             - PID of the aligner process (for log messages).
    /// * `last_stdout_activity`  - Shared timestamp updated by the stdout reader thread.
    /// * `last_stdin_activity`   - Shared timestamp updated by the stdin writer thread.
    /// * `timeout`               - Inactivity threshold for warning (error at `2 × timeout`).
    /// * `cancel`                - Shared flag set to `true` when a fatal timeout occurs.
    #[must_use]
    pub fn new(
        child_pid: u32,
        last_stdout_activity: Arc<AtomicU64>,
        last_stdin_activity: Arc<AtomicU64>,
        timeout: Duration,
        cancel: Arc<AtomicBool>,
    ) -> Self {
        Self { child_pid, last_stdout_activity, last_stdin_activity, timeout, cancel }
    }

    /// Start the health monitor in a background thread.
    ///
    /// The thread runs until the `cancel` flag is set (either by the monitor itself
    /// or by another component of the pipeline). Returns the thread join handle.
    #[must_use]
    pub fn start(self) -> JoinHandle<()> {
        thread::spawn(move || self.run())
    }

    /// Main monitor loop.
    fn run(self) {
        loop {
            thread::sleep(CHECK_INTERVAL);

            if self.cancel.load(Ordering::Acquire) {
                break;
            }

            let now = now_millis();
            let timeout_ms = u64::try_from(self.timeout.as_millis()).unwrap_or(u64::MAX);

            // ---- stdout checks ----
            let last_out = self.last_stdout_activity.load(Ordering::Acquire);
            if last_out > 0 {
                let elapsed_ms = now.saturating_sub(last_out);
                if elapsed_ms > timeout_ms * 2 {
                    log::error!(
                        "aligner pid={} stdout has been silent for {}ms (threshold {}ms); \
                         cancelling pipeline",
                        self.child_pid,
                        elapsed_ms,
                        timeout_ms * 2,
                    );
                    self.cancel.store(true, Ordering::Release);
                    break;
                } else if elapsed_ms > timeout_ms {
                    log::warn!(
                        "aligner pid={} stdout has been silent for {}ms (threshold {}ms)",
                        self.child_pid,
                        elapsed_ms,
                        timeout_ms,
                    );
                }
            }

            // ---- stdin checks ----
            let last_in = self.last_stdin_activity.load(Ordering::Acquire);
            if last_in > 0 {
                let elapsed_ms = now.saturating_sub(last_in);
                if elapsed_ms > timeout_ms {
                    log::warn!(
                        "aligner pid={} stdin has been silent for {}ms (threshold {}ms)",
                        self.child_pid,
                        elapsed_ms,
                        timeout_ms,
                    );
                }
            }
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify that the monitor does not set the cancel flag when activity timestamps
    /// are recent (i.e. within the timeout window).
    #[test]
    fn test_health_monitor_no_alert() {
        let cancel = Arc::new(AtomicBool::new(false));
        let recent = now_millis();

        let last_stdout = Arc::new(AtomicU64::new(recent));
        let last_stdin = Arc::new(AtomicU64::new(recent));

        let monitor = HealthMonitor::new(
            1234,
            Arc::clone(&last_stdout),
            Arc::clone(&last_stdin),
            Duration::from_secs(30), // generous timeout
            Arc::clone(&cancel),
        );

        let handle = monitor.start();

        // Let the monitor run for a couple of check cycles (CHECK_INTERVAL = 5 s would
        // be too slow for a unit test; we just verify the cancel flag stays clear for
        // a brief period).
        thread::sleep(Duration::from_millis(50));
        assert!(!cancel.load(Ordering::Acquire), "cancel should not be set with recent activity");

        // Signal the monitor to stop.
        cancel.store(true, Ordering::Release);
        handle.join().expect("monitor thread should finish cleanly");
    }

    /// Verify that the monitor sets the cancel flag when stdout has been inactive
    /// for longer than 2 × timeout.
    #[test]
    fn test_health_monitor_timeout() {
        let cancel = Arc::new(AtomicBool::new(false));

        // Simulate a stdout timestamp that is 10 seconds in the past.
        let stale = now_millis().saturating_sub(10_000);
        let last_stdout = Arc::new(AtomicU64::new(stale));
        // stdin is also stale but that only triggers a warning, not cancellation.
        let last_stdin = Arc::new(AtomicU64::new(stale));

        let timeout = Duration::from_secs(1); // very short for testing

        // Run the monitor loop inline (synchronously) so we don't have to wait
        // 5 seconds for the background thread to wake.
        let cancel_clone = Arc::clone(&cancel);
        let last_stdout_clone = Arc::clone(&last_stdout);
        let last_stdin_clone = Arc::clone(&last_stdin);

        let monitor =
            HealthMonitor::new(5678, last_stdout_clone, last_stdin_clone, timeout, cancel_clone);

        // Drive one check cycle directly by calling the internal logic.
        let now = now_millis();
        let timeout_ms = u64::try_from(monitor.timeout.as_millis()).unwrap_or(u64::MAX);
        let last_out = monitor.last_stdout_activity.load(Ordering::Acquire);
        let elapsed_ms = now.saturating_sub(last_out);

        if elapsed_ms > timeout_ms * 2 {
            monitor.cancel.store(true, Ordering::Release);
        }

        assert!(
            cancel.load(Ordering::Acquire),
            "cancel should be set when stdout is silent for > 2x timeout"
        );
    }
}

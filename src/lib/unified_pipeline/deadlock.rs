//! Deadlock detection and recovery for the unified pipeline.
//!
//! This module provides a watchdog mechanism that detects when the pipeline is stuck
//! (no progress for a configurable timeout) and logs detailed diagnostic info.
//! Optionally can recover by doubling queue memory limits.
//!
//! # Problem
//!
//! Despite non-blocking design with held-items pattern, deadlocks can still occur when:
//! - All threads have held items but can't push due to memory limits
//! - Reorder buffers fill up faster than they can drain
//! - Circular dependencies between backpressure conditions
//!
//! # Solution
//!
//! A progress tracking watchdog that:
//! 1. Tracks successful queue operations via per-queue timestamps
//! 2. Monitors for lack of progress over a timeout period (default: 10s)
//! 3. Default behavior: Logs detailed diagnostic info
//! 4. Optional recovery: With `--deadlock-recover`, uses adaptive limit management

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

/// Configuration for deadlock detection and recovery.
#[derive(Debug, Clone)]
pub struct DeadlockConfig {
    /// Timeout in seconds for deadlock detection (0 = disabled).
    pub timeout_secs: u64,
    /// Whether automatic recovery is enabled.
    pub recover_enabled: bool,
}

impl Default for DeadlockConfig {
    fn default() -> Self {
        Self {
            timeout_secs: 10, // Default 10 seconds
            recover_enabled: false,
        }
    }
}

impl DeadlockConfig {
    /// Create a new deadlock configuration.
    #[must_use]
    pub fn new(timeout_secs: u64, recover_enabled: bool) -> Self {
        Self { timeout_secs, recover_enabled }
    }

    /// Create a disabled configuration (no detection).
    #[must_use]
    pub fn disabled() -> Self {
        Self { timeout_secs: 0, recover_enabled: false }
    }

    /// Returns true if deadlock detection is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.timeout_secs > 0
    }
}

/// State for tracking pipeline progress and deadlock recovery.
///
/// This struct uses atomic operations for lock-free progress tracking from
/// multiple worker threads.
#[derive(Debug)]
pub struct DeadlockState {
    // Per-queue progress tracking (last activity timestamp in seconds since UNIX epoch)
    /// Read -> Q1 push timestamp
    pub q1_last_push: AtomicU64,
    /// Q1 -> Decompress pop timestamp
    pub q1_last_pop: AtomicU64,
    /// Decompress -> Q2 push timestamp
    pub q2_last_push: AtomicU64,
    /// Q2 -> `FindBoundaries` pop timestamp
    pub q2_last_pop: AtomicU64,
    /// `FindBoundaries` -> Q2b push timestamp (BAM pipeline)
    pub q2b_last_push: AtomicU64,
    /// Q2b -> Decode pop timestamp (BAM pipeline)
    pub q2b_last_pop: AtomicU64,
    /// `FindBoundaries` -> Q2.5 push timestamp (FASTQ pipeline)
    pub q2_5_last_push: AtomicU64,
    /// Q2.5 -> Parse pop timestamp (FASTQ pipeline)
    pub q2_5_last_pop: AtomicU64,
    /// Decode -> Q3 push timestamp
    pub q3_last_push: AtomicU64,
    /// Q3 -> Group pop timestamp
    pub q3_last_pop: AtomicU64,
    /// Group -> Q4 push timestamp
    pub q4_last_push: AtomicU64,
    /// Q4 -> Process pop timestamp
    pub q4_last_pop: AtomicU64,
    /// Process -> Q5 push timestamp
    pub q5_last_push: AtomicU64,
    /// Q5 -> Serialize pop timestamp
    pub q5_last_pop: AtomicU64,
    /// Serialize -> Q6 push timestamp
    pub q6_last_push: AtomicU64,
    /// Q6 -> Compress pop timestamp
    pub q6_last_pop: AtomicU64,
    /// Compress -> Q7 push timestamp
    pub q7_last_push: AtomicU64,
    /// Q7 -> Write pop timestamp
    pub q7_last_pop: AtomicU64,

    // Global progress tracking
    /// Unix timestamp of last progress on any queue.
    pub last_progress_time: AtomicU64,
    /// When last recovery happened (for restoration timer).
    pub last_recovery_time: AtomicU64,

    // Recovery state
    /// Current memory limit (modified by recovery). 0 = unlimited.
    pub current_memory_limit: AtomicU64,
    /// Original memory limit from config.
    pub original_memory_limit: u64,
    /// Lowest limit that worked (anti-oscillation).
    pub stable_memory_limit: AtomicU64,
    /// Did we deadlock at current limit?
    pub deadlock_at_current_limit: AtomicBool,

    // Configuration
    /// Deadlock detection timeout in seconds (0 = disabled).
    pub timeout_secs: u64,
    /// Whether automatic recovery is enabled.
    pub recover_enabled: bool,
}

impl DeadlockState {
    /// Create new deadlock state with the given configuration.
    #[must_use]
    pub fn new(config: &DeadlockConfig, original_memory_limit: u64) -> Self {
        let now = now_secs();
        Self {
            q1_last_push: AtomicU64::new(now),
            q1_last_pop: AtomicU64::new(now),
            q2_last_push: AtomicU64::new(now),
            q2_last_pop: AtomicU64::new(now),
            q2b_last_push: AtomicU64::new(now),
            q2b_last_pop: AtomicU64::new(now),
            q2_5_last_push: AtomicU64::new(now),
            q2_5_last_pop: AtomicU64::new(now),
            q3_last_push: AtomicU64::new(now),
            q3_last_pop: AtomicU64::new(now),
            q4_last_push: AtomicU64::new(now),
            q4_last_pop: AtomicU64::new(now),
            q5_last_push: AtomicU64::new(now),
            q5_last_pop: AtomicU64::new(now),
            q6_last_push: AtomicU64::new(now),
            q6_last_pop: AtomicU64::new(now),
            q7_last_push: AtomicU64::new(now),
            q7_last_pop: AtomicU64::new(now),
            last_progress_time: AtomicU64::new(now),
            last_recovery_time: AtomicU64::new(0),
            current_memory_limit: AtomicU64::new(original_memory_limit),
            original_memory_limit,
            stable_memory_limit: AtomicU64::new(0),
            deadlock_at_current_limit: AtomicBool::new(false),
            timeout_secs: config.timeout_secs,
            recover_enabled: config.recover_enabled,
        }
    }

    /// Create disabled state (no detection or tracking).
    #[must_use]
    pub fn disabled() -> Self {
        Self::new(&DeadlockConfig::disabled(), 0)
    }

    /// Returns true if deadlock detection is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.timeout_secs > 0
    }

    /// Get the current memory limit (for use in backpressure checks).
    #[must_use]
    pub fn get_memory_limit(&self) -> u64 {
        self.current_memory_limit.load(Ordering::Relaxed)
    }

    // ========== Progress Recording Methods ==========

    /// Record Q1 push progress.
    #[inline]
    pub fn record_q1_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q1_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q1 pop progress.
    #[inline]
    pub fn record_q1_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q1_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q2 push progress.
    #[inline]
    pub fn record_q2_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q2_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q2 pop progress.
    #[inline]
    pub fn record_q2_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q2_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q2b push progress.
    #[inline]
    pub fn record_q2b_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q2b_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q2b pop progress.
    #[inline]
    pub fn record_q2b_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q2b_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q2.5 push progress (FASTQ boundaries queue).
    #[inline]
    pub fn record_q2_5_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q2_5_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q2.5 pop progress (FASTQ boundaries queue).
    #[inline]
    pub fn record_q2_5_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q2_5_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q3 push progress.
    #[inline]
    pub fn record_q3_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q3_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q3 pop progress.
    #[inline]
    pub fn record_q3_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q3_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q4 push progress.
    #[inline]
    pub fn record_q4_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q4_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q4 pop progress.
    #[inline]
    pub fn record_q4_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q4_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q5 push progress.
    #[inline]
    pub fn record_q5_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q5_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q5 pop progress.
    #[inline]
    pub fn record_q5_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q5_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q6 push progress.
    #[inline]
    pub fn record_q6_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q6_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q6 pop progress.
    #[inline]
    pub fn record_q6_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q6_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q7 push progress.
    #[inline]
    pub fn record_q7_push(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q7_last_push.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }

    /// Record Q7 pop progress.
    #[inline]
    pub fn record_q7_pop(&self) {
        if self.is_enabled() {
            let now = now_secs();
            self.q7_last_pop.store(now, Ordering::Relaxed);
            self.last_progress_time.store(now, Ordering::Relaxed);
        }
    }
}

/// Get current time in seconds since UNIX epoch.
#[inline]
#[must_use]
pub fn now_secs() -> u64 {
    SystemTime::now().duration_since(UNIX_EPOCH).map_or(0, |d| d.as_secs())
}

/// Snapshot of queue state for diagnostics.
#[derive(Debug, Clone)]
pub struct QueueSnapshot {
    pub q1_len: usize,
    pub q2_len: usize,
    pub q2b_len: usize,
    pub q3_len: usize,
    pub q4_len: usize,
    pub q5_len: usize,
    pub q6_len: usize,
    pub q7_len: usize,
    pub q2_reorder_mem: u64,
    pub q3_reorder_mem: u64,
    pub memory_limit: u64,
    pub read_done: bool,
    pub group_done: bool,
    pub draining: bool,
}

/// Check for deadlock and attempt recovery if configured.
///
/// This function should be called periodically from the monitor thread (every ~1 second).
///
/// # Returns
///
/// - `DeadlockAction::None` if no deadlock detected
/// - `DeadlockAction::Detected` if deadlock detected but recovery disabled
/// - `DeadlockAction::Recovered(new_limit)` if recovery was performed
pub fn check_deadlock_and_restore(
    deadlock_state: &DeadlockState,
    snapshot: &QueueSnapshot,
) -> DeadlockAction {
    if !deadlock_state.is_enabled() {
        return DeadlockAction::None;
    }

    let now = now_secs();
    let last_progress = deadlock_state.last_progress_time.load(Ordering::Relaxed);
    let timeout = deadlock_state.timeout_secs;

    // Check for deadlock (no progress for timeout period)
    if now.saturating_sub(last_progress) >= timeout {
        return handle_deadlock(deadlock_state, snapshot, now);
    }

    // Check for restoration opportunity (30s of sustained progress after recovery)
    if deadlock_state.recover_enabled {
        try_restore_limits(deadlock_state, now);
    }

    DeadlockAction::None
}

/// Result of deadlock check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeadlockAction {
    /// No deadlock detected.
    None,
    /// Deadlock detected but recovery disabled.
    Detected,
    /// Deadlock detected and recovery performed with new limit.
    Recovered(u64),
}

/// Handle a detected deadlock.
fn handle_deadlock(
    deadlock_state: &DeadlockState,
    snapshot: &QueueSnapshot,
    now: u64,
) -> DeadlockAction {
    let current = deadlock_state.current_memory_limit.load(Ordering::Relaxed);

    // Log diagnostics
    log::warn!("DEADLOCK DETECTED: No progress for {}s", deadlock_state.timeout_secs);
    log_queue_state(deadlock_state, snapshot);

    let action = if deadlock_state.recover_enabled {
        // Mark that we deadlocked at this limit (for anti-oscillation)
        deadlock_state.deadlock_at_current_limit.store(true, Ordering::Relaxed);

        // Update stable limit - current limit caused deadlock, so stable is 2x current
        if current > 0 {
            let new_stable = current.saturating_mul(2);
            let old_stable = deadlock_state.stable_memory_limit.load(Ordering::Relaxed);
            if new_stable > old_stable {
                deadlock_state.stable_memory_limit.store(new_stable, Ordering::Relaxed);
                log::warn!("  Stable limit updated: {} -> {} bytes", old_stable, new_stable);
            }
        }

        // Double the limit
        let new_limit = if current == 0 {
            0 // Already unbounded
        } else {
            current.saturating_mul(2)
        };

        let final_limit = if new_limit == 0
            || new_limit > 8u64.saturating_mul(deadlock_state.original_memory_limit)
        {
            // Cap at 8x or unbind
            log::warn!("  Recovery: Unbinding limits (unlimited)");
            0
        } else {
            log::warn!("  Recovery: {} -> {} bytes (2x)", current, new_limit);
            new_limit
        };

        deadlock_state.current_memory_limit.store(final_limit, Ordering::SeqCst);
        deadlock_state.last_recovery_time.store(now, Ordering::Relaxed);

        DeadlockAction::Recovered(final_limit)
    } else {
        log::warn!("  Recovery disabled - use --deadlock-recover to enable automatic recovery");
        DeadlockAction::Detected
    };

    // Reset progress timer to avoid repeated detections
    deadlock_state.last_progress_time.store(now, Ordering::Relaxed);

    action
}

/// Try to restore limits after sustained progress.
fn try_restore_limits(deadlock_state: &DeadlockState, now: u64) {
    let last_recovery = deadlock_state.last_recovery_time.load(Ordering::Relaxed);

    // Only restore if we've had 30s of progress since last recovery
    if last_recovery == 0 || now.saturating_sub(last_recovery) < 30 {
        return;
    }

    let current = deadlock_state.current_memory_limit.load(Ordering::Relaxed);
    let original = deadlock_state.original_memory_limit;
    let stable = deadlock_state.stable_memory_limit.load(Ordering::Relaxed);

    if current == 0 {
        // Unbounded - try to restore to 8x original or stable, whichever is higher
        let target = (original.saturating_mul(8)).max(stable);
        if target > 0 {
            log::info!("  Restoring: unlimited -> {} bytes", target);
            deadlock_state.current_memory_limit.store(target, Ordering::SeqCst);
            deadlock_state.last_recovery_time.store(now, Ordering::Relaxed);
            deadlock_state.deadlock_at_current_limit.store(false, Ordering::Relaxed);
        }
        return;
    }

    // Only restore if current limit is above original and above stable limit
    if current <= original || (stable > 0 && current <= stable) {
        return; // Already at minimum safe limit
    }

    // Halve the limit (but not below original or stable)
    let new_limit = (current / 2).max(original).max(stable);
    if new_limit < current {
        log::info!("  Restoring: {} -> {} bytes (halving)", current, new_limit);
        deadlock_state.current_memory_limit.store(new_limit, Ordering::SeqCst);
        deadlock_state.last_recovery_time.store(now, Ordering::Relaxed);
        deadlock_state.deadlock_at_current_limit.store(false, Ordering::Relaxed);
    }
}

/// Log detailed queue state for diagnostics.
fn log_queue_state(deadlock_state: &DeadlockState, snapshot: &QueueSnapshot) {
    let now = now_secs();

    // Queue depths
    log::warn!(
        "  Queue depths: Q1={} Q2={} Q2b={} Q3={} Q4={} Q5={} Q6={} Q7={}",
        snapshot.q1_len,
        snapshot.q2_len,
        snapshot.q2b_len,
        snapshot.q3_len,
        snapshot.q4_len,
        snapshot.q5_len,
        snapshot.q6_len,
        snapshot.q7_len,
    );

    // Per-queue stall detection (seconds since last activity)
    log::warn!(
        "  Q1: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q1_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q1_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q2: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q2_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q2_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q2b: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q2b_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q2b_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q2.5: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q2_5_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q2_5_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q3: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q3_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q3_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q4: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q4_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q4_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q5: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q5_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q5_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q6: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q6_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q6_last_pop.load(Ordering::Relaxed))
    );
    log::warn!(
        "  Q7: push={}s ago, pop={}s ago",
        now.saturating_sub(deadlock_state.q7_last_push.load(Ordering::Relaxed)),
        now.saturating_sub(deadlock_state.q7_last_pop.load(Ordering::Relaxed))
    );

    // Memory state
    log::warn!(
        "  Memory: Q2 reorder={:.1}MB, Q3 reorder={:.1}MB, limit={:.1}MB",
        snapshot.q2_reorder_mem as f64 / 1_000_000.0,
        snapshot.q3_reorder_mem as f64 / 1_000_000.0,
        snapshot.memory_limit as f64 / 1_000_000.0,
    );

    // Pipeline state flags
    log::warn!(
        "  State: read_done={}, group_done={}, draining={}",
        snapshot.read_done,
        snapshot.group_done,
        snapshot.draining,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    // ========== Progress Tracking Tests ==========

    #[test]
    fn test_progress_timestamp_updates_on_push() {
        let config = DeadlockConfig::new(10, false);
        let state = DeadlockState::new(&config, 512 * 1024 * 1024);

        let initial = state.last_progress_time.load(Ordering::Relaxed);

        // Sleep briefly to ensure time advances
        std::thread::sleep(std::time::Duration::from_millis(10));

        state.record_q1_push();
        let after_push = state.last_progress_time.load(Ordering::Relaxed);

        // Progress time should be updated (or at least not decreased)
        assert!(after_push >= initial);
    }

    #[test]
    fn test_progress_timestamp_updates_on_pop() {
        let config = DeadlockConfig::new(10, false);
        let state = DeadlockState::new(&config, 512 * 1024 * 1024);

        let initial = state.last_progress_time.load(Ordering::Relaxed);

        std::thread::sleep(std::time::Duration::from_millis(10));

        state.record_q1_pop();
        let after_pop = state.last_progress_time.load(Ordering::Relaxed);

        assert!(after_pop >= initial);
    }

    #[test]
    fn test_per_queue_timestamps_independent() {
        let config = DeadlockConfig::new(10, false);
        let state = DeadlockState::new(&config, 512 * 1024 * 1024);

        // Record on different queues
        state.record_q1_push();
        let q1_time = state.q1_last_push.load(Ordering::Relaxed);

        std::thread::sleep(std::time::Duration::from_millis(10));

        state.record_q3_push();
        let q3_time = state.q3_last_push.load(Ordering::Relaxed);

        // Q1 timestamp should not have changed
        assert_eq!(state.q1_last_push.load(Ordering::Relaxed), q1_time);
        // Q3 timestamp should be updated
        assert!(q3_time >= q1_time);
    }

    #[test]
    fn test_disabled_state_skips_tracking() {
        let state = DeadlockState::disabled();

        let initial = state.last_progress_time.load(Ordering::Relaxed);
        state.record_q1_push();
        let after = state.last_progress_time.load(Ordering::Relaxed);

        // When disabled, timestamps don't update
        assert_eq!(initial, after);
    }

    // ========== Deadlock Detection Tests ==========

    #[test]
    fn test_deadlock_detected_after_timeout() {
        let config = DeadlockConfig::new(1, false); // 1 second timeout
        let state = DeadlockState::new(&config, 512 * 1024 * 1024);

        // Set last_progress_time to 5 seconds ago
        let now = now_secs();
        state.last_progress_time.store(now.saturating_sub(5), Ordering::Relaxed);

        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 10,
            q2b_len: 5,
            q3_len: 10,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 100_000_000,
            q3_reorder_mem: 200_000_000,
            memory_limit: 512_000_000,
            read_done: false,
            group_done: false,
            draining: false,
        };

        let action = check_deadlock_and_restore(&state, &snapshot);
        assert_eq!(action, DeadlockAction::Detected);
    }

    #[test]
    fn test_no_detection_before_timeout() {
        let config = DeadlockConfig::new(10, false); // 10 second timeout
        let state = DeadlockState::new(&config, 512 * 1024 * 1024);

        // Progress happened just now (state initialized with current time)
        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 0,
            q2b_len: 0,
            q3_len: 0,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 0,
            q3_reorder_mem: 0,
            memory_limit: 0,
            read_done: false,
            group_done: false,
            draining: false,
        };

        let action = check_deadlock_and_restore(&state, &snapshot);
        assert_eq!(action, DeadlockAction::None);
    }

    #[test]
    fn test_detection_disabled_returns_none() {
        let state = DeadlockState::disabled();

        // Set an old progress time
        let now = now_secs();
        state.last_progress_time.store(now.saturating_sub(100), Ordering::Relaxed);

        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 0,
            q2b_len: 0,
            q3_len: 0,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 0,
            q3_reorder_mem: 0,
            memory_limit: 0,
            read_done: false,
            group_done: false,
            draining: false,
        };

        let action = check_deadlock_and_restore(&state, &snapshot);
        assert_eq!(action, DeadlockAction::None);
    }

    // ========== Recovery Tests ==========

    #[test]
    fn test_recovery_doubles_limit() {
        let config = DeadlockConfig::new(1, true); // recovery enabled
        let original_limit = 512 * 1024 * 1024; // 512 MB
        let state = DeadlockState::new(&config, original_limit);

        // Simulate deadlock
        let now = now_secs();
        state.last_progress_time.store(now.saturating_sub(5), Ordering::Relaxed);

        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 10,
            q2b_len: 5,
            q3_len: 10,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 400_000_000,
            q3_reorder_mem: 500_000_000,
            memory_limit: original_limit,
            read_done: false,
            group_done: false,
            draining: false,
        };

        let action = check_deadlock_and_restore(&state, &snapshot);

        // Should have doubled the limit
        if let DeadlockAction::Recovered(new_limit) = action {
            assert_eq!(new_limit, original_limit * 2);
        } else {
            panic!("Expected Recovered action, got {:?}", action);
        }
    }

    #[test]
    fn test_recovery_caps_at_8x() {
        let config = DeadlockConfig::new(1, true);
        let original_limit = 512 * 1024 * 1024; // 512 MB
        let state = DeadlockState::new(&config, original_limit);

        // Set current limit to 4GB (8x 512MB)
        state.current_memory_limit.store(original_limit * 8, Ordering::Relaxed);

        let now = now_secs();
        state.last_progress_time.store(now.saturating_sub(5), Ordering::Relaxed);

        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 10,
            q2b_len: 5,
            q3_len: 10,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 0,
            q3_reorder_mem: 0,
            memory_limit: original_limit * 8,
            read_done: false,
            group_done: false,
            draining: false,
        };

        let action = check_deadlock_and_restore(&state, &snapshot);

        // Should unbind (0 = unlimited)
        if let DeadlockAction::Recovered(new_limit) = action {
            assert_eq!(new_limit, 0);
        } else {
            panic!("Expected Recovered action, got {:?}", action);
        }
    }

    #[test]
    fn test_recovery_updates_stable_limit() {
        let config = DeadlockConfig::new(1, true);
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        // Initial stable limit should be 0
        assert_eq!(state.stable_memory_limit.load(Ordering::Relaxed), 0);

        let now = now_secs();
        state.last_progress_time.store(now.saturating_sub(5), Ordering::Relaxed);

        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 0,
            q2b_len: 0,
            q3_len: 0,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 0,
            q3_reorder_mem: 0,
            memory_limit: original_limit,
            read_done: false,
            group_done: false,
            draining: false,
        };

        check_deadlock_and_restore(&state, &snapshot);

        // Stable limit should be 2x the failed limit
        let stable = state.stable_memory_limit.load(Ordering::Relaxed);
        assert_eq!(stable, original_limit * 2);
    }

    #[test]
    fn test_recovery_disabled_only_logs() {
        let config = DeadlockConfig::new(1, false); // recovery disabled
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        let now = now_secs();
        state.last_progress_time.store(now.saturating_sub(5), Ordering::Relaxed);

        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 0,
            q2b_len: 0,
            q3_len: 0,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 0,
            q3_reorder_mem: 0,
            memory_limit: original_limit,
            read_done: false,
            group_done: false,
            draining: false,
        };

        let action = check_deadlock_and_restore(&state, &snapshot);

        // Should detect but not recover
        assert_eq!(action, DeadlockAction::Detected);
        // Limit should be unchanged
        assert_eq!(state.current_memory_limit.load(Ordering::Relaxed), original_limit);
    }

    // ========== Restoration Tests ==========

    #[test]
    fn test_restoration_after_30s_progress() {
        let config = DeadlockConfig::new(10, true);
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        // Set current limit to 2x original (simulating previous recovery)
        state.current_memory_limit.store(original_limit * 2, Ordering::SeqCst);

        // Set recovery time to 35s ago
        let now = now_secs();
        state.last_recovery_time.store(now.saturating_sub(35), Ordering::Relaxed);

        // Call try_restore_limits directly
        try_restore_limits(&state, now);

        // Should have halved the limit back to original
        let new_limit = state.current_memory_limit.load(Ordering::Relaxed);
        assert_eq!(new_limit, original_limit);
    }

    #[test]
    fn test_no_restoration_before_30s() {
        let config = DeadlockConfig::new(10, true);
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        // Set current limit to 2x original
        let elevated_limit = original_limit * 2;
        state.current_memory_limit.store(elevated_limit, Ordering::SeqCst);

        // Set recovery time to 20s ago (not enough time)
        let now = now_secs();
        state.last_recovery_time.store(now.saturating_sub(20), Ordering::Relaxed);

        try_restore_limits(&state, now);

        // Limit should be unchanged
        assert_eq!(state.current_memory_limit.load(Ordering::Relaxed), elevated_limit);
    }

    #[test]
    fn test_restoration_respects_stable_limit() {
        let config = DeadlockConfig::new(10, true);
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        // Set stable limit to 1GB (higher than original)
        let stable_limit = 1024 * 1024 * 1024;
        state.stable_memory_limit.store(stable_limit, Ordering::Relaxed);

        // Set current limit to 2GB
        let current_limit = 2 * 1024 * 1024 * 1024;
        state.current_memory_limit.store(current_limit, Ordering::SeqCst);

        // Set recovery time to 35s ago
        let now = now_secs();
        state.last_recovery_time.store(now.saturating_sub(35), Ordering::Relaxed);

        try_restore_limits(&state, now);

        // Should have halved to 1GB, but not below stable limit
        let new_limit = state.current_memory_limit.load(Ordering::Relaxed);
        assert_eq!(new_limit, stable_limit);
    }

    #[test]
    fn test_restoration_from_unbounded() {
        let config = DeadlockConfig::new(10, true);
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        // Set current limit to 0 (unbounded)
        state.current_memory_limit.store(0, Ordering::SeqCst);

        // Set recovery time to 35s ago
        let now = now_secs();
        state.last_recovery_time.store(now.saturating_sub(35), Ordering::Relaxed);

        try_restore_limits(&state, now);

        // Should restore to 8x original
        let new_limit = state.current_memory_limit.load(Ordering::Relaxed);
        assert_eq!(new_limit, original_limit * 8);
    }

    // ========== Anti-Oscillation Tests ==========

    #[test]
    fn test_stable_limit_never_decreases() {
        let config = DeadlockConfig::new(1, true);
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        // Set stable limit to 1GB
        let stable_limit = 1024 * 1024 * 1024;
        state.stable_memory_limit.store(stable_limit, Ordering::Relaxed);

        // Try to update with a lower value (shouldn't work via deadlock detection)
        let now = now_secs();
        state.last_progress_time.store(now.saturating_sub(5), Ordering::Relaxed);

        // Set current limit to something that would result in lower stable
        state.current_memory_limit.store(256 * 1024 * 1024, Ordering::Relaxed);

        let snapshot = QueueSnapshot {
            q1_len: 0,
            q2_len: 0,
            q2b_len: 0,
            q3_len: 0,
            q4_len: 0,
            q5_len: 0,
            q6_len: 0,
            q7_len: 0,
            q2_reorder_mem: 0,
            q3_reorder_mem: 0,
            memory_limit: 256 * 1024 * 1024,
            read_done: false,
            group_done: false,
            draining: false,
        };

        check_deadlock_and_restore(&state, &snapshot);

        // Stable limit should still be 1GB (never decreased)
        let current_stable = state.stable_memory_limit.load(Ordering::Relaxed);
        assert_eq!(current_stable, stable_limit);
    }

    // ========== Configuration Tests ==========

    #[test]
    fn test_config_default() {
        let config = DeadlockConfig::default();
        assert_eq!(config.timeout_secs, 10);
        assert!(!config.recover_enabled);
        assert!(config.is_enabled());
    }

    #[test]
    fn test_config_disabled() {
        let config = DeadlockConfig::disabled();
        assert_eq!(config.timeout_secs, 0);
        assert!(!config.recover_enabled);
        assert!(!config.is_enabled());
    }

    #[test]
    fn test_get_memory_limit() {
        let config = DeadlockConfig::new(10, true);
        let original_limit = 512 * 1024 * 1024;
        let state = DeadlockState::new(&config, original_limit);

        assert_eq!(state.get_memory_limit(), original_limit);

        // Modify and verify
        state.current_memory_limit.store(1024 * 1024 * 1024, Ordering::Relaxed);
        assert_eq!(state.get_memory_limit(), 1024 * 1024 * 1024);
    }
}

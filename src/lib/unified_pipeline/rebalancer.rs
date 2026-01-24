//! Dynamic memory rebalancer for pipeline queues.
//!
//! This module provides functionality to dynamically adjust memory allocations
//! across pipeline queues based on observed demand, preventing starvation and
//! optimizing throughput.
//!
//! # Algorithm
//!
//! The rebalancer collects statistics every epoch (default: 1 second) and
//! adjusts allocations based on:
//! - Peak memory usage (indicates demand)
//! - Time spent blocked (indicates starvation)
//! - Utilization ratio (high utilization suggests need for more)
//!
//! Changes are applied gradually (max 20% per epoch) to avoid oscillation.

use super::queue::QueueStats;
use std::time::Duration;

/// Configuration for the dynamic rebalancer.
#[derive(Debug, Clone)]
pub struct RebalancerConfig {
    /// Total memory budget across all queues.
    pub total_budget: u64,
    /// How often to rebalance (default: 1 second).
    pub epoch_duration: Duration,
    /// Maximum change per epoch as a fraction (0.2 = 20%).
    pub max_change_fraction: f64,
    /// Minimum allocation per queue (bytes).
    pub min_per_queue: u64,
}

impl Default for RebalancerConfig {
    fn default() -> Self {
        Self {
            total_budget: 4 * 1024 * 1024 * 1024, // 4GB
            epoch_duration: Duration::from_secs(1),
            max_change_fraction: 0.2,
            min_per_queue: 64 * 1024 * 1024, // 64MB
        }
    }
}

/// Manages dynamic memory allocation across queues.
///
/// The rebalancer monitors queue statistics and adjusts allocations to
/// optimize throughput while staying within the total memory budget.
///
/// # Example
///
/// ```ignore
/// let config = RebalancerConfig {
///     total_budget: 2 * 1024 * 1024 * 1024, // 2GB
///     ..Default::default()
/// };
/// let names = vec!["q2_reorder", "q3_reorder", "q7_reorder"];
/// let initial = vec![500_000_000, 1_000_000_000, 500_000_000];
/// let mut rebalancer = DynamicRebalancer::new(config, names, initial);
///
/// // After collecting stats
/// let stats = vec![queue1.collect_stats(), queue2.collect_stats(), queue3.collect_stats()];
/// let new_limits = rebalancer.rebalance(&stats);
/// ```
pub struct DynamicRebalancer {
    config: RebalancerConfig,
    current_allocations: Vec<u64>,
    queue_names: Vec<&'static str>,
}

impl DynamicRebalancer {
    /// Create a new rebalancer with the given configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - Rebalancer configuration
    /// * `queue_names` - Names of queues (for logging)
    /// * `initial_allocations` - Starting allocations for each queue
    ///
    /// # Panics
    ///
    /// Panics if `queue_names` and `initial_allocations` have different lengths.
    #[must_use]
    pub fn new(
        config: RebalancerConfig,
        queue_names: Vec<&'static str>,
        initial_allocations: Vec<u64>,
    ) -> Self {
        assert_eq!(queue_names.len(), initial_allocations.len());
        Self { config, current_allocations: initial_allocations, queue_names }
    }

    /// Compute new allocations based on observed stats.
    ///
    /// # Algorithm
    ///
    /// 1. Compute demand for each queue based on peak usage and starvation
    /// 2. Allocate proportionally to demand
    /// 3. Apply gradual change limits to avoid oscillation
    /// 4. Enforce minimum allocation per queue
    /// 5. Normalize to exactly match total budget
    ///
    /// # Arguments
    ///
    /// * `stats` - Statistics collected from each queue for this epoch
    ///
    /// # Returns
    ///
    /// New memory limits for each queue (in same order as `queue_names`).
    pub fn rebalance(&mut self, stats: &[QueueStats]) -> Vec<u64> {
        assert_eq!(stats.len(), self.current_allocations.len());

        // Compute demand for each queue
        let demands: Vec<f64> = stats
            .iter()
            .zip(self.current_allocations.iter())
            .map(|(s, &current_alloc)| {
                // Base demand is peak usage (or minimum)
                let base_demand = s.peak_bytes.max(self.config.min_per_queue) as f64;

                // Starvation factor: if blocked, increase demand
                // Scale by blocked time in seconds
                let starvation_factor = 1.0 + (s.time_blocked_ms as f64 / 1000.0);

                // Utilization factor: if using most of allocation, likely needs more
                let utilization = s.peak_bytes as f64 / current_alloc.max(1) as f64;
                let utilization_factor = if utilization > 0.9 { 1.5 } else { 1.0 };

                base_demand * starvation_factor * utilization_factor
            })
            .collect();

        let total_demand: f64 = demands.iter().sum();

        // Compute ideal allocations proportional to demand
        let ideal: Vec<u64> = demands
            .iter()
            .map(|&d| {
                let fraction = d / total_demand.max(1.0);
                (fraction * self.config.total_budget as f64) as u64
            })
            .collect();

        // Apply gradual change limit and minimum
        let new_allocations: Vec<u64> = self
            .current_allocations
            .iter()
            .zip(ideal.iter())
            .map(|(&current, &target)| {
                let max_change = (current as f64 * self.config.max_change_fraction) as u64;
                let max_change = max_change.max(self.config.min_per_queue / 4);

                let new = if target > current {
                    current.saturating_add(max_change.min(target - current))
                } else {
                    current.saturating_sub(max_change.min(current - target))
                };

                new.max(self.config.min_per_queue)
            })
            .collect();

        // Normalize to exactly match budget
        let sum: u64 = new_allocations.iter().sum();
        let normalized: Vec<u64> = if sum > 0 {
            new_allocations
                .iter()
                .map(|&a| ((a as f64 / sum as f64) * self.config.total_budget as f64) as u64)
                .collect()
        } else {
            // Fallback: equal distribution
            let per_queue = self.config.total_budget / new_allocations.len() as u64;
            vec![per_queue; new_allocations.len()]
        };

        self.current_allocations.clone_from(&normalized);
        normalized
    }

    /// Get current allocations.
    #[must_use]
    pub fn current_allocations(&self) -> &[u64] {
        &self.current_allocations
    }

    /// Get queue names.
    #[must_use]
    pub fn queue_names(&self) -> &[&'static str] {
        &self.queue_names
    }

    /// Log current state (for debugging/monitoring).
    pub fn log_state(&self, stats: &[QueueStats]) {
        log::debug!("Queue memory rebalance:");
        for (i, name) in self.queue_names.iter().enumerate() {
            let alloc = self.current_allocations[i];
            let stat = &stats[i];
            log::debug!(
                "  {}: {}/{} MB (peak: {} MB, blocked: {}ms)",
                name,
                stat.peak_bytes / (1024 * 1024),
                alloc / (1024 * 1024),
                stat.peak_bytes / (1024 * 1024),
                stat.time_blocked_ms
            );
        }
    }
}

/// Initial memory allocation across queues.
#[derive(Debug, Clone)]
pub struct InitialAllocation {
    /// Allocation for Q2 reorder buffer (Decompress → `FindBoundaries`).
    pub q2_reorder: u64,
    /// Allocation for Q3 reorder buffer (Decode → Group).
    pub q3_reorder: u64,
    /// Allocation for Q7 reorder buffer (Compress → Write).
    pub q7_reorder: u64,
    /// Allocation for array queues (FIFO queues between steps).
    pub array_queues: u64,
}

/// Get initial memory allocations based on command and strategy.
///
/// These percentages are empirically derived from benchmarks to provide
/// good starting points. The dynamic rebalancer will adjust from here.
///
/// # Arguments
///
/// * `command` - The fgumi command (e.g., "group", "simplex", "duplex")
/// * `strategy` - Optional strategy for the command (e.g., "identity", "edit")
/// * `total_budget` - Total memory budget in bytes
///
/// # Returns
///
/// Initial allocations for each queue type.
#[must_use]
pub fn initial_allocation_for_command(
    command: &str,
    strategy: Option<&str>,
    total_budget: u64,
) -> InitialAllocation {
    // Percentages derived from benchmark training
    // Format: (q2_pct, q3_pct, q7_pct, other_pct) - must sum to 100
    let (q2_pct, q3_pct, q7_pct, other_pct) = match (command, strategy) {
        // Group strategies - Q3 typically dominates due to grouping accumulation
        ("group", Some("identity")) => (20, 50, 15, 15),
        ("group", Some("edit")) => (25, 45, 15, 15),
        ("group", Some("adjacency")) => (30, 40, 15, 15),
        ("group", Some("paired")) => (25, 45, 15, 15),
        ("group", _) => (25, 45, 15, 15), // Default for group

        // Consensus commands
        ("simplex", _) => (25, 45, 15, 15),
        ("duplex", _) => (25, 45, 15, 15),
        ("codec", _) => (25, 45, 15, 15),

        // Filter - may have larger output
        ("filter", _) => (25, 40, 20, 15),

        // Correct - similar to group
        ("correct", _) => (30, 40, 15, 15),

        // Extract (FASTQ pipeline) - different queue structure
        ("extract", _) => (35, 35, 15, 15),

        // Default
        _ => (25, 40, 15, 20),
    };

    InitialAllocation {
        q2_reorder: total_budget * q2_pct / 100,
        q3_reorder: total_budget * q3_pct / 100,
        q7_reorder: total_budget * q7_pct / 100,
        array_queues: total_budget * other_pct / 100,
    }
}

/// Parse human-friendly memory values like "4GB", "512MB", "1024M".
///
/// Supports the following suffixes (case-insensitive):
/// - GB, G: Gigabytes
/// - MB, M: Megabytes
/// - KB, K: Kilobytes
/// - No suffix: Bytes
///
/// # Arguments
///
/// * `s` - The string to parse (e.g., "4GB", "512MB", "1024")
///
/// # Returns
///
/// The value in bytes, or an error message.
///
/// # Errors
///
/// Returns an error if:
/// - The value cannot be parsed as a number
/// - The value is below the minimum (256MB)
pub fn parse_memory_limit(s: &str) -> Result<u64, String> {
    let s = s.trim().to_uppercase();

    let (num_str, multiplier) = if s.ends_with("GB") {
        (s.trim_end_matches("GB"), 1024 * 1024 * 1024)
    } else if s.ends_with('G') {
        (s.trim_end_matches('G'), 1024 * 1024 * 1024)
    } else if s.ends_with("MB") {
        (s.trim_end_matches("MB"), 1024 * 1024)
    } else if s.ends_with('M') {
        (s.trim_end_matches('M'), 1024 * 1024)
    } else if s.ends_with("KB") {
        (s.trim_end_matches("KB"), 1024)
    } else if s.ends_with('K') {
        (s.trim_end_matches('K'), 1024)
    } else {
        // Assume bytes
        (s.as_str(), 1)
    };

    let num: f64 = num_str.trim().parse().map_err(|_| format!("Invalid memory value: {}", s))?;

    let bytes = (num * f64::from(multiplier)) as u64;

    const MIN_BYTES: u64 = 256 * 1024 * 1024; // 256MB
    if bytes < MIN_BYTES {
        return Err(format!(
            "Queue memory limit {}MB is too low (minimum: 256MB)",
            bytes / (1024 * 1024)
        ));
    }

    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_memory_limit() {
        assert_eq!(parse_memory_limit("4GB").unwrap(), 4 * 1024 * 1024 * 1024);
        assert_eq!(parse_memory_limit("4G").unwrap(), 4 * 1024 * 1024 * 1024);
        assert_eq!(parse_memory_limit("512MB").unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_memory_limit("512M").unwrap(), 512 * 1024 * 1024);
        assert_eq!(parse_memory_limit("1gb").unwrap(), 1024 * 1024 * 1024);
        assert_eq!(parse_memory_limit("  2GB  ").unwrap(), 2 * 1024 * 1024 * 1024);
        assert_eq!(parse_memory_limit("1.5GB").unwrap(), (1.5 * 1024.0 * 1024.0 * 1024.0) as u64);
    }

    #[test]
    fn test_parse_memory_limit_minimum() {
        assert!(parse_memory_limit("100MB").is_err());
        assert!(parse_memory_limit("256MB").is_ok());
    }

    #[test]
    fn test_initial_allocation() {
        let alloc = initial_allocation_for_command("group", Some("adjacency"), 1_000_000_000);
        assert_eq!(alloc.q2_reorder, 300_000_000); // 30%
        assert_eq!(alloc.q3_reorder, 400_000_000); // 40%
        assert_eq!(alloc.q7_reorder, 150_000_000); // 15%
        assert_eq!(alloc.array_queues, 150_000_000); // 15%
    }

    #[test]
    fn test_rebalancer_basic() {
        let config = RebalancerConfig {
            total_budget: 1_000_000_000, // 1GB
            min_per_queue: 100_000_000,  // 100MB
            ..Default::default()
        };
        let names = vec!["q1", "q2", "q3"];
        let initial = vec![333_333_333, 333_333_333, 333_333_334];
        let mut rebalancer = DynamicRebalancer::new(config, names, initial);

        // Simulate q2 being starved
        let stats = vec![
            QueueStats { avg_bytes: 100_000_000, peak_bytes: 200_000_000, time_blocked_ms: 0 },
            QueueStats {
                avg_bytes: 300_000_000,
                peak_bytes: 330_000_000,
                time_blocked_ms: 500, // Blocked for 500ms
            },
            QueueStats { avg_bytes: 100_000_000, peak_bytes: 150_000_000, time_blocked_ms: 0 },
        ];

        let new_allocs = rebalancer.rebalance(&stats);

        // q2 should get more allocation due to starvation
        assert!(new_allocs[1] > new_allocs[0]);
        assert!(new_allocs[1] > new_allocs[2]);

        // Total should be close to budget (allow for rounding)
        let total: u64 = new_allocs.iter().sum();
        let budget: u64 = 1_000_000_000;
        assert!(
            total >= budget - 10 && total <= budget + 10,
            "Total {} should be close to budget {}",
            total,
            budget
        );
    }

    #[test]
    fn test_rebalancer_gradual_change() {
        let config = RebalancerConfig {
            total_budget: 1_000_000_000,
            max_change_fraction: 0.2,
            min_per_queue: 100_000_000,
            ..Default::default()
        };
        let names = vec!["q1", "q2"];
        let initial = vec![500_000_000, 500_000_000];
        let mut rebalancer = DynamicRebalancer::new(config, names, initial);

        // Simulate q1 needing much more
        let stats = vec![
            QueueStats {
                avg_bytes: 450_000_000,
                peak_bytes: 500_000_000,
                time_blocked_ms: 1000, // Very starved
            },
            QueueStats { avg_bytes: 50_000_000, peak_bytes: 100_000_000, time_blocked_ms: 0 },
        ];

        let new_allocs = rebalancer.rebalance(&stats);

        // Change should be limited to ~20%
        let change = (new_allocs[0] as i64 - 500_000_000i64).abs();
        assert!(change <= 150_000_000); // Allow some slack for normalization
    }
}

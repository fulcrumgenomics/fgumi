//! Cgroup-aware system resource detection.
//!
//! Provides [`detect_total_memory`] and [`detect_cpu_count`] which correctly
//! report resource limits when running inside Docker or Kubernetes containers.
//! Both functions take the minimum of the cgroup limit and the host value so
//! they behave correctly on bare-metal too.
//!
//! # Memory (`detect_total_memory`)
//!
//! Uses `sysinfo::System::cgroup_limits()` which reads:
//! - cgroup v2: `/sys/fs/cgroup/memory.max`
//! - cgroup v1: `/sys/fs/cgroup/memory/memory.limit_in_bytes`
//!
//! Falls back to `System::total_memory()` (physical RAM) on macOS and when
//! no cgroup limit is configured.
//!
//! # CPU (`detect_cpu_count`)
//!
//! Uses `std::thread::available_parallelism`, which reads the CFS quota:
//! - cgroup v2: `cpu.max`
//! - cgroup v1: `cpu.cfs_quota_us` / `cpu.cfs_period_us`
//!
//! Falls back to the platform core count when no quota is configured.
//! Note: CFS quota is the *hard* limit set by `--cpus`; it is not the same as
//! `cpu.shares` which is a soft scheduling weight.

/// Returns the effective total memory available to this process in bytes.
///
/// Checks cgroup memory limits for container environments, falling back to
/// physical RAM on bare-metal or macOS. Always returns `min(cgroup, physical)`.
#[must_use]
pub fn detect_total_memory() -> usize {
    let mut system = sysinfo::System::new();
    system.refresh_memory();
    let physical = system.total_memory();
    let bytes = system.cgroup_limits().map_or(physical, |c| c.total_memory.min(physical));
    // Saturate at usize::MAX / 2 rather than usize::MAX on 32-bit platforms so the
    // downstream `budget > total` overflow check in `resolve_memory_budget` can fire
    // correctly (no value can exceed usize::MAX, so using it as the fallback renders
    // the check dead).
    usize::try_from(bytes).unwrap_or(usize::MAX / 2)
}

/// Returns the number of logical CPUs available to this process.
///
/// Respects cgroup CPU quotas (set by `--cpus` in Docker or resource limits in
/// Kubernetes), falling back to the platform core count. Returns at least 1.
#[must_use]
pub fn detect_cpu_count() -> usize {
    std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_total_memory_nonzero() {
        let total = detect_total_memory();
        assert!(total > 0, "expected non-zero total memory, got {total}");
    }

    #[test]
    fn test_detect_total_memory_bounded_by_physical() {
        let total = detect_total_memory();
        let mut system = sysinfo::System::new();
        system.refresh_memory();
        let physical = usize::try_from(system.total_memory()).unwrap_or(usize::MAX);
        assert!(total <= physical, "cgroup-limited total {total} exceeded physical {physical}");
    }

    #[test]
    fn test_detect_cpu_count_at_least_one() {
        assert!(detect_cpu_count() >= 1);
    }

    #[test]
    fn test_detect_cpu_count_reasonable() {
        // No machine has more than 65536 logical CPUs (yet).
        assert!(detect_cpu_count() <= 65536);
    }
}

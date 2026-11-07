//! A monotonic clock for telemetry that does not read the clock-source hardware.
//!
//! [`std::time::Instant::now`] reads `CLOCK_MONOTONIC`. The vDSO normally serves
//! that from the TSC in a few tens of nanoseconds, but when the kernel's clock
//! source cannot be read from user space (`hpet` or `acpi_pm`, which a kernel
//! falls back to after marking the TSC unstable, as on some VMs) every read is
//! a system call that touches the device, costing microseconds. Timing every
//! step dispatch with it then turns idle polling into system time.
//!
//! On Linux, [`coarse_now_ns`] reads `CLOCK_MONOTONIC_COARSE` instead. The
//! vDSO serves that from the time the kernel stored at its last tick, whatever
//! the clock source, so it never leaves user space. Its resolution is one
//! kernel tick (1 to 10 ms, depending on `CONFIG_HZ`): a single short interval
//! reads as zero or as one tick. Summed over many intervals the rounding
//! averages out, which is what busy/idle totals need, but a per-interval value
//! is not meaningful below the tick. On other platforms `Instant` is cheap and
//! [`coarse_now_ns`] reads it, so it has `Instant`'s resolution there.

#[cfg(not(target_os = "linux"))]
use std::sync::OnceLock;
#[cfg(not(target_os = "linux"))]
use std::time::Instant;

/// The current time on the coarse monotonic clock, in nanoseconds from an
/// arbitrary, process-wide origin. Only differences between two values are
/// meaningful.
#[cfg(target_os = "linux")]
#[inline]
#[must_use]
pub(crate) fn coarse_now_ns() -> u64 {
    let ts = rustix::time::clock_gettime(rustix::time::ClockId::MonotonicCoarse);
    // A monotonic clock never reports a negative time.
    let secs = u64::try_from(ts.tv_sec).unwrap_or(0);
    let nanos = u64::try_from(ts.tv_nsec).unwrap_or(0);
    secs.saturating_mul(1_000_000_000).saturating_add(nanos)
}

/// The current time on the coarse monotonic clock, in nanoseconds from an
/// arbitrary, process-wide origin. Only differences between two values are
/// meaningful.
#[cfg(not(target_os = "linux"))]
#[inline]
#[must_use]
pub(crate) fn coarse_now_ns() -> u64 {
    static ORIGIN: OnceLock<Instant> = OnceLock::new();
    let origin = *ORIGIN.get_or_init(Instant::now);
    u64::try_from(origin.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use super::coarse_now_ns;

    #[test]
    fn now_never_goes_backwards() {
        let mut prev = coarse_now_ns();
        for _ in 0..10_000 {
            let next = coarse_now_ns();
            assert!(next >= prev, "{next} < {prev}");
            prev = next;
        }
    }

    #[test]
    fn elapsed_tracks_a_sleep_within_one_tick() {
        // A tick is at most 10 ms (`CONFIG_HZ=100`), so a 50 ms sleep reads as
        // at least 40 ms on the coarse clock.
        let start = coarse_now_ns();
        std::thread::sleep(Duration::from_millis(50));
        let elapsed = coarse_now_ns().saturating_sub(start);
        assert!(elapsed >= 40_000_000, "elapsed {elapsed} ns after a 50 ms sleep");
    }
}

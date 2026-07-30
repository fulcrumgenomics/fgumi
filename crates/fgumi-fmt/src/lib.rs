#![deny(unsafe_code)]

//! Human-readable formatting helpers shared across the fgumi workspace.
//!
//! These are needed by crates at several layers that do not — and should not — depend on one
//! another, so they live in a leaf crate near the bottom of the graph. Anything that only
//! turns a number into a string for humans to read belongs here rather than in a crate that
//! would drag `clap`, `sysinfo`, or BAM I/O along with it.

use num_format::{Locale, ToFormattedString};
use std::time::Duration;

/// Format an integer with comma thousands-separators (e.g. `1234567` → `"1,234,567"`).
///
/// # Examples
///
/// ```
/// use fgumi_fmt::format_count;
///
/// assert_eq!(format_count(1_234_567), "1,234,567");
/// assert_eq!(format_count(123), "123");
/// assert_eq!(format_count(0), "0");
/// ```
#[must_use]
pub fn format_count(n: u64) -> String {
    n.to_formatted_string(&Locale::en)
}

/// Formats a duration in human-readable form (e.g. `"45s"`, `"2m 15s"`, `"1h 30m"`).
///
/// Renders at most two units (whole seconds below a minute, minutes-and-seconds below an
/// hour, hours-and-minutes above); the finer unit is dropped when it is zero.
///
/// # Examples
///
/// ```
/// use fgumi_fmt::format_duration;
/// use std::time::Duration;
///
/// assert_eq!(format_duration(Duration::from_secs(45)), "45s");
/// assert_eq!(format_duration(Duration::from_secs(135)), "2m 15s");
/// assert_eq!(format_duration(Duration::from_secs(5400)), "1h 30m");
/// ```
#[must_use]
pub fn format_duration(duration: Duration) -> String {
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

/// Formats a throughput as items per second, falling back to items per minute below one
/// item per second.
///
/// Durations under a millisecond are reported as `count` items/s rather than dividing by a
/// near-zero elapsed time, which would produce a meaningless rate.
///
/// # Examples
///
/// ```
/// use fgumi_fmt::format_rate;
/// use std::time::Duration;
///
/// assert_eq!(format_rate(600, Duration::from_secs(60)), "10 items/s");
/// // Below one item per second the units switch to items/min.
/// assert_eq!(format_rate(30, Duration::from_secs(60)), "30.0 items/min");
/// ```
#[must_use]
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn format_rate(count: u64, duration: Duration) -> String {
    let secs = duration.as_secs_f64();
    if secs < 0.001 {
        return format!("{} items/s", format_count(count));
    }

    let rate = count as f64 / secs;
    if rate >= 1.0 {
        format!("{} items/s", format_count(rate as u64))
    } else {
        let items_per_min = count as f64 / (secs / 60.0);
        format!("{items_per_min:.1} items/min")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::zero(0, "0")]
    #[case::single(1, "1")]
    #[case::three_digits(123, "123")]
    #[case::first_comma(1_234, "1,234")]
    #[case::two_groups(1_234_567, "1,234,567")]
    #[case::round_million(1_000_000, "1,000,000")]
    #[case::max(u64::MAX, "18,446,744,073,709,551,615")]
    fn format_count_inserts_commas(#[case] n: u64, #[case] expected: &str) {
        assert_eq!(format_count(n), expected);
    }

    #[rstest]
    #[case::zero(0, "0s")]
    #[case::sub_minute(45, "45s")]
    #[case::exact_minute(60, "1m")]
    #[case::minutes_seconds(135, "2m 15s")]
    #[case::exact_hour(3600, "1h")]
    #[case::hours_minutes(5400, "1h 30m")]
    fn format_duration_human_readable(#[case] secs: u64, #[case] expected: &str) {
        assert_eq!(format_duration(Duration::from_secs(secs)), expected);
    }

    /// `format_rate` has three regimes: a sub-millisecond guard that avoids dividing by
    /// ~zero, a normal items/s path, and an items/min path below one item per second.
    #[rstest]
    #[case::sub_millisecond_guard(5, Duration::from_micros(100), "5 items/s")]
    #[case::exactly_one_per_second(1, Duration::from_secs(1), "1 items/s")]
    #[case::thousands_per_second(1000, Duration::from_secs(1), "1,000 items/s")]
    #[case::below_one_per_second(30, Duration::from_secs(60), "30.0 items/min")]
    #[case::far_below_one_per_second(1, Duration::from_secs(120), "0.5 items/min")]
    fn format_rate_switches_units_below_one_per_second(
        #[case] count: u64,
        #[case] duration: Duration,
        #[case] expected: &str,
    ) {
        assert_eq!(format_rate(count, duration), expected);
    }
}

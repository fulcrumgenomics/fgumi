#![deny(unsafe_code)]

//! Small, dependency-free human-readable formatting helpers shared across the fgumi
//! workspace.
//!
//! `format_count` and `format_duration` were each copied into more than one crate
//! (`fgumi-cli-common`, `fgumi-metrics`, `fgumi-bam-io`). Those crates sit at different
//! layers and do not depend on one another, so this dependency-free leaf crate is the
//! shared home; each consumer re-exports or imports from here instead of carrying its own
//! copy.

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
    let s = n.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    let num_commas = if len > 3 { (len - 1) / 3 } else { 0 };
    let mut result = String::with_capacity(len + num_commas);
    for (i, &byte) in bytes.iter().enumerate() {
        // A comma precedes every group of three trailing digits, but never leads the
        // string: `len - i` is the number of digits remaining, so a comma falls before
        // position `i` exactly when that count is a positive multiple of three.
        if i > 0 && (len - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(byte as char);
    }
    result
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
}

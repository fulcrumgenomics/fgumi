//! When the UMI-matching index is used instead of a linear scan.

use std::fmt;
use std::str::FromStr;

/// Default minimum number of unique UMIs to use N-gram/BK-tree index for candidate search.
/// Below this threshold, linear scan is used. Based on benchmarks: index provides ~10%
/// speedup even at 100 UMIs/position with negligible construction overhead.
pub const DEFAULT_INDEX_THRESHOLD: usize = 100;

/// Value of `--index-threshold`: when to match UMIs through the N-gram/BK-tree index
/// rather than by comparing every pair.
///
/// The index is a pure optimization — it narrows the candidate set, it does not change
/// which UMIs group together — so every variant yields the same grouping and differs only
/// in how much work reaching it takes.
///
/// [`MinUmis`](Self::MinUmis) is a *minimum group size*, not a switch, and that has caught
/// people out: `MinUmis(0)` admits every group, so a plain `0` turns the index fully **on**
/// — the opposite of what the CLI help claimed for a long time. [`Always`](Self::Always) and
/// [`Never`](Self::Never) let that intent be stated in a word rather than encoded in a
/// number no one can read back. `Never` is not otherwise reachable: suppressing the index
/// through `MinUmis` alone would take a value larger than any position group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexThreshold {
    /// Index every position group, whatever its size.
    Always,
    /// Never index; always compare all UMI pairs.
    Never,
    /// Index a position group once it holds at least this many distinct UMIs.
    MinUmis(usize),
}

impl IndexThreshold {
    /// Whether a position group of `num_umis` distinct UMIs is matched through the index.
    #[must_use]
    pub fn admits(self, num_umis: usize) -> bool {
        match self {
            Self::Always => true,
            Self::Never => false,
            Self::MinUmis(min) => num_umis >= min,
        }
    }

    /// Whether this value *asserts* that indexing must happen, rather than merely tuning
    /// when it starts.
    ///
    /// Only the `always` keyword asserts. A number — including `0`, which admits every
    /// group exactly as `Always` does — is a tuning knob allowed to end up inert, which is
    /// what lets the default `100` coexist with strategies that never index. Callers use
    /// this to reject a request they cannot honour instead of silently ignoring it.
    #[must_use]
    pub fn demands_indexing(self) -> bool {
        matches!(self, Self::Always)
    }
}

impl Default for IndexThreshold {
    fn default() -> Self {
        Self::MinUmis(DEFAULT_INDEX_THRESHOLD)
    }
}

impl From<usize> for IndexThreshold {
    fn from(min_umis: usize) -> Self {
        Self::MinUmis(min_umis)
    }
}

impl FromStr for IndexThreshold {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("always") {
            Ok(Self::Always)
        } else if s.eq_ignore_ascii_case("never") {
            Ok(Self::Never)
        } else {
            s.parse::<usize>().map(Self::MinUmis).map_err(|e| {
                format!("expected 'always', 'never', or a non-negative integer, got '{s}': {e}")
            })
        }
    }
}

impl fmt::Display for IndexThreshold {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Always => f.write_str("always"),
            Self::Never => f.write_str("never"),
            Self::MinUmis(min) => write!(f, "{min}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// `--index-threshold` accepts the two keywords and any non-negative integer.
    #[rstest]
    #[case::always("always", IndexThreshold::Always)]
    #[case::always_mixed_case("Always", IndexThreshold::Always)]
    #[case::never("never", IndexThreshold::Never)]
    #[case::never_upper("NEVER", IndexThreshold::Never)]
    #[case::zero("0", IndexThreshold::MinUmis(0))]
    #[case::default_value("100", IndexThreshold::MinUmis(100))]
    fn test_index_threshold_parses(#[case] input: &str, #[case] expected: IndexThreshold) {
        assert_eq!(input.parse::<IndexThreshold>(), Ok(expected));
    }

    /// Rejected input names the accepted forms, so the error is actionable.
    #[rstest]
    #[case::negative("-1")]
    #[case::word("off")]
    #[case::empty("")]
    #[case::float("1.5")]
    fn test_index_threshold_rejects_invalid(#[case] input: &str) {
        let err = input.parse::<IndexThreshold>().expect_err("must reject");
        assert!(
            err.contains("always") && err.contains("never"),
            "error should name the accepted forms, got: {err}"
        );
    }

    /// `admits` is the gate the assigners consult. `MinUmis` is a MINIMUM, so `0`
    /// admits every group — the reading the CLI help had backwards for a long time.
    #[rstest]
    #[case::always_admits_empty(IndexThreshold::Always, 0, true)]
    #[case::always_admits_huge(IndexThreshold::Always, 1_000_000, true)]
    #[case::never_rejects_huge(IndexThreshold::Never, 1_000_000, false)]
    #[case::zero_admits_singleton(IndexThreshold::MinUmis(0), 1, true)]
    #[case::below_minimum(IndexThreshold::MinUmis(100), 99, false)]
    #[case::at_minimum(IndexThreshold::MinUmis(100), 100, true)]
    #[case::above_minimum(IndexThreshold::MinUmis(100), 101, true)]
    fn test_index_threshold_admits(
        #[case] threshold: IndexThreshold,
        #[case] num_umis: usize,
        #[case] expected: bool,
    ) {
        assert_eq!(threshold.admits(num_umis), expected);
    }

    /// Only the `always` keyword asserts that indexing must happen. A number is a
    /// tuning knob allowed to end up inert, which is what lets the default `100`
    /// coexist with strategies that never index.
    #[rstest]
    #[case::always_demands(IndexThreshold::Always, true)]
    #[case::never_does_not(IndexThreshold::Never, false)]
    #[case::zero_does_not(IndexThreshold::MinUmis(0), false)]
    #[case::default_does_not(IndexThreshold::MinUmis(100), false)]
    fn test_only_always_demands_indexing(
        #[case] threshold: IndexThreshold,
        #[case] expected: bool,
    ) {
        assert_eq!(threshold.demands_indexing(), expected);
    }

    /// Values round-trip through `Display`, so logs and help text stay parseable.
    #[rstest]
    #[case::always(IndexThreshold::Always, "always")]
    #[case::never(IndexThreshold::Never, "never")]
    #[case::min_umis(IndexThreshold::MinUmis(100), "100")]
    fn test_index_threshold_displays_and_round_trips(
        #[case] threshold: IndexThreshold,
        #[case] expected: &str,
    ) {
        assert_eq!(threshold.to_string(), expected);
        assert_eq!(expected.parse::<IndexThreshold>(), Ok(threshold));
    }

    /// The default matches the documented `--index-threshold` default.
    #[test]
    fn test_default_is_the_documented_threshold() {
        assert_eq!(IndexThreshold::default(), IndexThreshold::MinUmis(DEFAULT_INDEX_THRESHOLD));
    }
}

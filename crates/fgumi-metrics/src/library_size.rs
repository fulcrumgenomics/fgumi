//! Lander–Waterman library-size estimation, a faithful port of Picard/HTSJDK
//! `DuplicationMetrics.estimateLibrarySize`. Pure and IO-free.

/// Estimate the number of unique molecules in a library from the observed
/// read-pair counts, using the Lander–Waterman / Picard model.
///
/// `n_pairs` is the total number of (observed) read pairs; `n_unique` is the
/// number of unique (non-duplicate) read pairs. Returns `None` when the
/// estimate is undefined: no pairs, no duplicates (`n_unique == n_pairs`), or
/// `n_unique >= n_pairs`.
///
/// This mirrors Picard's implementation bit-for-bit: it bisects a ratio `r`
/// (library size `= n_unique * r`) over `f(x) = c/x - 1 + exp(-n/x)` for 40
/// iterations after bracketing, and returns `n_unique * r` truncated to an
/// integer.
#[must_use]
#[expect(
    clippy::cast_precision_loss,
    reason = "read-pair counts stay well under 2^52 for any realistic library"
)]
#[expect(
    clippy::many_single_char_names,
    reason = "faithful port of Picard's estimateLibrarySize; names mirror the math (n, c, f, m, r, u)"
)]
pub fn estimate_library_size(n_pairs: u64, n_unique: u64) -> Option<u64> {
    // Undefined cases (Picard returns null / guards these).
    if n_pairs == 0 || n_unique == 0 || n_unique >= n_pairs {
        return None;
    }
    let n = n_pairs as f64;
    let c = n_unique as f64;

    // f(x) = c/x - 1 + exp(-n/x); x = c * ratio.
    let f = |x: f64| c / x - 1.0 + (-n / x).exp();

    let mut m = 1.0_f64;
    let mut big_m = 100.0_f64;

    // Picard bails if f(m*c) < 0 at the low end; treat as undefined.
    if f(m * c) < 0.0 {
        return None;
    }
    // Expand the upper bracket until f(M*c) < 0.
    while f(big_m * c) >= 0.0 {
        big_m *= 10.0;
    }
    // 40-step bisection over the ratio.
    for _ in 0..40 {
        let r = f64::midpoint(m, big_m);
        let u = f(r * c);
        if u == 0.0 {
            break;
        } else if u > 0.0 {
            m = r;
        } else {
            big_m = r;
        }
    }
    let estimate = c * f64::midpoint(m, big_m);
    #[expect(
        clippy::cast_possible_truncation,
        reason = "Picard truncates the final ratio*c product to an integer library size"
    )]
    #[expect(clippy::cast_sign_loss, reason = "estimate is always non-negative by construction")]
    let estimate = estimate as u64;
    Some(estimate)
}

#[cfg(test)]
mod tests {
    use super::estimate_library_size;
    use rstest::rstest;

    #[rstest]
    #[case::typical_1m(1_000_000, 900_000, Some(4_660_793))]
    #[case::typical_1k(1_000, 900, Some(4_660))]
    #[case::half_unique(500_000, 250_000, Some(313_750))]
    #[case::high_complexity(10_000, 9_500, Some(96_638))]
    #[case::no_duplicates(100, 100, None)] // unique == pairs
    #[case::unique_exceeds(100, 200, None)] // unique > pairs
    #[case::zero_pairs(0, 0, None)]
    #[case::zero_unique(1000, 0, None)]
    fn matches_picard_reference(
        #[case] n_pairs: u64,
        #[case] n_unique: u64,
        #[case] expected: Option<u64>,
    ) {
        let got = estimate_library_size(n_pairs, n_unique);
        match (got, expected) {
            (Some(g), Some(e)) => assert!(
                (i128::from(g) - i128::from(e)).abs() <= 1,
                "n_pairs={n_pairs} n_unique={n_unique}: got {g}, expected {e}"
            ),
            (a, b) => assert_eq!(a, b, "n_pairs={n_pairs} n_unique={n_unique}"),
        }
    }
}

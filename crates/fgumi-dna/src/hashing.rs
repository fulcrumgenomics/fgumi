//! Deterministic, contention-free hashing for hot per-item maps.

use ahash::RandomState;

/// A fixed-seed [`RandomState`] for `AHashMap`s created in hot per-item loops
/// (per-position-group, per-assign call).
///
/// `AHashMap::new()` / `::default()` derives its state from [`RandomState::new`], which bumps a
/// process-global atomic counter on every construction. Under UMI grouping's
/// many-small-maps-per-thread pattern that counter's cache line ping-pongs across cores and
/// dominates CPU; giving each map a fixed seed skips the counter entirely (measured: `group`
/// throughput +12-37% at 16 threads on a 53M-record BAM, every output byte-identical).
///
/// The seed is arbitrary and reproducible across runs. Output never depends on it: every consumer
/// sorts before assigning molecule IDs or mints IDs in input order, so map iteration order cannot
/// reach the result. This is a throughput fix, not a determinism fix.
///
/// This is the single home for the fixed seed table; `fgumi_sort::cb_hasher` and the per-crate
/// `hashing` re-exports all delegate here.
#[inline]
#[must_use]
pub fn deterministic_state() -> RandomState {
    // Arbitrary fixed seeds — chosen for uniqueness, not cryptographic strength.
    RandomState::with_seeds(
        0xa1b2_c3d4_e5f6_0718,
        0x9182_7364_5546_3728,
        0xfede_dcba_0987_6543,
        0x0011_2233_4455_6677,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two independently-constructed `deterministic_state()` values must hash an identical key to
    /// the identical value. This is the fixed-seed contract the throughput fix relies on: it holds
    /// for `RandomState::with_seeds` but would fail if `deterministic_state` were ever swapped back
    /// to `RandomState::new()`, whose seeds differ on every construction.
    #[rstest::rstest]
    #[case::empty(b"")]
    #[case::short(b"ACGT")]
    #[case::umi(b"ACGT-TGCA")]
    #[case::long(b"the quick brown fox jumps over the lazy dog")]
    fn deterministic_state_hashes_key_identically_across_instances(#[case] key: &[u8]) {
        assert_eq!(deterministic_state().hash_one(key), deterministic_state().hash_one(key));
    }
}

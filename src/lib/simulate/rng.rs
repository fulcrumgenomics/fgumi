//! Random number generator utilities for simulation.
//!
//! Provides seeded RNG creation for reproducible simulations.

use rand::SeedableRng;
use rand::rngs::StdRng;

/// Create a random number generator, optionally seeded for reproducibility.
///
/// # Arguments
///
/// * `seed` - Optional seed value. If `None`, uses OS entropy for randomness.
///
/// # Returns
///
/// A `StdRng` instance that can be used for all simulation randomness.
///
/// # Examples
///
/// ```
/// use fgumi_lib::simulate::create_rng;
///
/// // Reproducible simulation
/// let mut rng1 = create_rng(Some(42));
/// let mut rng2 = create_rng(Some(42));
/// // rng1 and rng2 will produce identical sequences
///
/// // Random simulation (different each run)
/// let mut rng3 = create_rng(None);
/// ```
#[must_use]
pub fn create_rng(seed: Option<u64>) -> StdRng {
    match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => rand::make_rng(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::RngExt;

    #[test]
    fn test_seeded_rng_reproducible() {
        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(42));

        let values1: Vec<u64> = (0..10).map(|_| rng1.random()).collect();
        let values2: Vec<u64> = (0..10).map(|_| rng2.random()).collect();

        assert_eq!(values1, values2);
    }

    #[test]
    fn test_different_seeds_different_values() {
        let mut rng1 = create_rng(Some(42));
        let mut rng2 = create_rng(Some(43));

        let values1: Vec<u64> = (0..10).map(|_| rng1.random()).collect();
        let values2: Vec<u64> = (0..10).map(|_| rng2.random()).collect();

        assert_ne!(values1, values2);
    }

    #[test]
    fn test_unseeded_rng_works() {
        let mut rng = create_rng(None);
        let _value: u64 = rng.random();
        // Just verify it doesn't panic
    }
}

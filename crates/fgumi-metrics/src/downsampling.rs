//! Deterministic downsampling via Murmur3 hashing.
//!
//! Provides a fixed set of downsampling fractions and a hash function that maps
//! read names to a value in [0, 1], enabling deterministic subsetting of reads.

use murmur3::murmur3_32;

/// Standard downsampling fractions: 5%, 10%, 15%, ..., 100%.
pub const DOWNSAMPLING_FRACTIONS: [f64; 20] = [
    0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,
    0.85, 0.90, 0.95, 1.00,
];

/// Computes a hash value normalized to the [0, 1] range using Murmur3.
///
/// Hashes the read name using the Murmur3 32-bit algorithm with seed 42, then
/// normalizes the result to a floating-point value between 0 and 1. This is used
/// for deterministic downsampling where reads are assigned to fractions based on
/// their hash value. Matches the Scala implementation's hashing approach.
///
/// # Arguments
///
/// * `read_name` - The read name string to hash
///
/// # Returns
///
/// A floating-point value in the range [0, 1].
#[must_use]
pub fn compute_hash_fraction(read_name: &str) -> f64 {
    // Use Murmur3 with seed 42 (matching Scala implementation)
    let hash = murmur3_32(&mut std::io::Cursor::new(read_name.as_bytes()), 42).unwrap_or(0);

    // Scala implementation uses Int (i32), takes absolute value, then normalizes
    // Convert u32 to i32 first to match Scala's behavior
    // unsigned_abs(i32::MIN) == 2147483648 > i32::MAX, so clamp to avoid values > 1.0
    let positive_hash = f64::from(hash.cast_signed().unsigned_abs());
    (positive_hash / f64::from(i32::MAX)).min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash_fraction_in_range() {
        for name in ["read1", "read2", "AAACCC", "some:long:read:name:123"] {
            let frac = compute_hash_fraction(name);
            assert!(
                (0.0..=1.0).contains(&frac),
                "hash fraction {frac} for '{name}' is out of [0, 1]"
            );
        }
    }

    #[test]
    fn test_hash_fraction_deterministic() {
        let a = compute_hash_fraction("read1");
        let b = compute_hash_fraction("read1");
        assert!((a - b).abs() < f64::EPSILON, "hash fraction should be deterministic");
    }

    #[test]
    fn test_hash_fraction_differs_for_different_names() {
        let a = compute_hash_fraction("read1");
        let b = compute_hash_fraction("read2");
        assert!(
            (a - b).abs() > f64::EPSILON,
            "different names should (almost certainly) produce different fractions"
        );
    }

    #[test]
    fn test_downsampling_fractions_count() {
        assert_eq!(DOWNSAMPLING_FRACTIONS.len(), 20);
        assert!((DOWNSAMPLING_FRACTIONS[0] - 0.05).abs() < f64::EPSILON);
        assert!((DOWNSAMPLING_FRACTIONS[19] - 1.00).abs() < f64::EPSILON);
    }
}

//! Typed errors for the `fgumi-umi` crate.
//!
//! Publicly-returned errors use [`enum@Error`] so that callers can pattern-match on
//! specific failure kinds. Internal helpers may still use `anyhow::Error` when
//! the added churn of typed errors is not worth it.

use thiserror::Error;

/// Errors returned by the public API of `fgumi-umi`.
#[derive(Debug, Error)]
pub enum Error {
    /// A paired UMI string did not have the expected `A-B` form (exactly one
    /// `-` separator).
    #[error(
        "UMI '{umi}' is not a valid paired UMI (expected format: 'A-B'{})",
        if *.multiple_separators { ", found multiple '-'" } else { "" }
    )]
    InvalidPairedUmi {
        /// The offending UMI string.
        umi: String,
        /// True if the UMI contained more than one `-` separator; false if it
        /// contained none.
        multiple_separators: bool,
    },
}

/// Convenience `Result` alias for [`enum@Error`].
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_variants_are_matchable() {
        let err = Error::InvalidPairedUmi { umi: "ACGT".to_string(), multiple_separators: false };
        assert!(matches!(err, Error::InvalidPairedUmi { multiple_separators: false, .. }));
    }

    #[test]
    fn display_formats_both_cases() {
        let missing =
            Error::InvalidPairedUmi { umi: "ACGT".to_string(), multiple_separators: false };
        assert!(missing.to_string().contains("ACGT"));
        assert!(!missing.to_string().contains("multiple"));

        let multiple =
            Error::InvalidPairedUmi { umi: "A-B-C".to_string(), multiple_separators: true };
        assert!(multiple.to_string().contains("A-B-C"));
        assert!(multiple.to_string().contains("multiple"));
    }
}

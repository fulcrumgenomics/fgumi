//! Typed errors for the `fgumi-consensus` crate.
//!
//! Publicly-returned errors use [`enum@Error`] so that callers can pattern-match on
//! specific failure kinds. Internal helpers may still use `anyhow::Error` when
//! the added churn of typed errors is not worth it.

use thiserror::Error;

/// Errors returned by the public API of `fgumi-consensus`.
#[derive(Debug, Error)]
pub enum Error {
    /// A configuration parameter was invalid. Covers bad `min_reads` vector,
    /// mis-sized UMI tag names, and other config-space validation failures.
    #[error("invalid consensus config: {0}")]
    InvalidConfig(String),

    /// A raw BAM record was malformed: too short, missing a required alignment
    /// position, out-of-bounds index access during consensus construction, etc.
    #[error("invalid BAM record: {0}")]
    InvalidRecord(String),

    /// A required SAM tag was missing from a record.
    #[error("missing required tag '{tag}' for read '{read_name}'")]
    MissingTag {
        /// The two-letter tag identifier (e.g. `"MI"`, `"RX"`).
        tag: String,
        /// The read name (best-effort, may be lossy UTF-8).
        read_name: String,
    },

    /// A SAM tag had an invalid value: non-UTF-8 bytes, malformed MI suffix, etc.
    #[error("invalid tag value: {0}")]
    InvalidTag(String),

    /// Consensus calling was attempted on an empty input.
    #[error("cannot create consensus from empty source reads")]
    EmptySourceReads,

    /// Codec duplex disagreement exceeded the configured threshold.
    #[error("codec duplex disagreement exceeded threshold: {0}")]
    DuplexDisagreement(String),

    /// Failure serializing a record back to raw BAM bytes. Captures the
    /// underlying `std::io::Error`.
    #[error("failed to encode BAM record: {source}")]
    Encode {
        /// The underlying I/O error from the encoder.
        #[from]
        source: std::io::Error,
    },
}

/// Convenience `Result` alias for [`enum@Error`].
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_variants_are_matchable() {
        let err = Error::InvalidConfig("bad threshold".into());
        assert!(matches!(err, Error::InvalidConfig(_)));

        let err = Error::EmptySourceReads;
        assert!(matches!(err, Error::EmptySourceReads));

        let err = Error::MissingTag { tag: "MI".into(), read_name: "r1".into() };
        assert!(matches!(err, Error::MissingTag { .. }));

        let err = Error::Encode { source: std::io::Error::other("boom") };
        assert!(matches!(err, Error::Encode { .. }));
    }

    #[test]
    fn io_error_auto_converts_via_from() {
        fn inner() -> Result<()> {
            Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "bad"))?;
            Ok(())
        }
        assert!(matches!(inner(), Err(Error::Encode { .. })));
    }
}
